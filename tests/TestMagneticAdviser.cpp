#include "RandomUtils.h"
#include <source_location>
#include "support/Painter.h"
#include "advisers/MagneticAdviser.h"
#include "processors/CircuitSimulatorInterface.h"
#include "processors/Inputs.h"
#include "TestingUtils.h"
#include "processors/Sweeper.h"
#include "physical_models/Impedance.h"
#include "converter_models/IsolatedBuck.h"
#include "Definitions.h"
#include <magic_enum.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

namespace {
    bool plot = false;
    std::string file_path = std::source_location::current().file_name();

    TEST_CASE("Test_MagneticAdviser_Filter_Available_Cores", "[adviser][magnetic-adviser]") {
        SKIP("Test needs investigation");
        clear_databases();
        settings.set_use_concentric_cores(false);
        settings.set_use_toroidal_cores(true);
        settings.set_use_only_cores_in_stock(false);

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {24, 78, 76};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 507026;
        double magnetizingInductance = 100e-6;
        double temperature = 25;
        WaveformLabel waveShape = WaveformLabel::TRIANGULAR;
        double peakToPeak = 100;
        double dutyCycle = 0.5;
        double dcCurrent = 0;

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        ImpedancePoint impedancePoint;
        impedancePoint.set_magnitude(5000);
        ImpedanceAtFrequency impedanceAtFrequency;
        impedanceAtFrequency.set_frequency(1e5);
        impedanceAtFrequency.set_impedance(impedancePoint);
        inputs.get_mutable_design_requirements().set_minimum_impedance(std::vector<ImpedanceAtFrequency>{impedanceAtFrequency});
        inputs.get_mutable_design_requirements().get_mutable_magnetizing_inductance().set_minimum(magnetizingInductance);
        inputs.get_mutable_design_requirements().get_mutable_magnetizing_inductance().set_nominal(std::nullopt);
        inputs.get_mutable_design_requirements().get_mutable_magnetizing_inductance().set_maximum(std::nullopt);

        OpenMagnetics::Mas masMagnetic;
        inputs.process();

        MagneticAdviser magneticAdviser;
        magneticAdviser.set_application(MAS::Application::INTERFERENCE_SUPPRESSION);
        magneticAdviser.set_core_mode(CoreAdviser::CoreAdviserModes::AVAILABLE_CORES);
        auto masMagnetics = magneticAdviser.get_advised_magnetic(inputs, 1);
        REQUIRE(masMagnetics.size() > 0);

        for (auto [masMagnetic, scoring] : masMagnetics) {
            auto impedance = Impedance().calculate_impedance(masMagnetic.get_mutable_magnetic(), 1e6);
            REQUIRE(abs(impedance) >= 5000);

            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            if (plot) {
                auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
                // {
                //     auto outFile = outputFilePath;
                //     std::string filename = "Test_MagneticAdviser_Filter_" + std::to_string(scoring) + ".svg";
                //     outFile.append(filename);
                //     settings.set_painter_include_fringing(false);
                //     Painter painter(outFile, true);

                //     painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
                //     painter.paint_core(masMagnetic.get_mutable_magnetic());
                //     painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
                //     painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
                //     painter.export_svg();
                // }
                {
                    auto impedanceSweep = Sweeper().sweep_impedance_over_frequency(masMagnetic.get_mutable_magnetic(), 10000, 400e6, 1000);

                    auto outFile = outputFilePath;

                    outFile.append("Test_MagneticAdviser_Filter_Sweep_Impedance_Over_Frequency_" + std::to_string(scoring) + ".svg");
                    std::filesystem::remove(outFile);
                    Painter painter(outFile, false, true);
                    #ifdef ENABLE_MATPLOTPP
                    painter.paint_curve(impedanceSweep);
                    #else
                        INFO("matplotplusplus disabled — skipping AdvancedPainter call");
                    #endif

                    #ifdef ENABLE_MATPLOTPP
                    painter.export_svg();
                    #else
                        INFO("matplotplusplus disabled — skipping AdvancedPainter call");
                    #endif

                }
            }
        }
        settings.reset();
    }

    TEST_CASE("Test_MagneticAdviser_Filter_Standard_Cores", "[adviser][magnetic-adviser]") {
        SKIP("Test needs investigation");
        clear_databases();
        settings.set_use_concentric_cores(false);
        settings.set_use_toroidal_cores(true);
        settings.set_use_only_cores_in_stock(false);

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {24, 78, 76};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 507026;
        double magnetizingInductance = 100e-6;
        double temperature = 25;
        WaveformLabel waveShape = WaveformLabel::TRIANGULAR;
        double peakToPeak = 100;
        double dutyCycle = 0.5;
        double dcCurrent = 0;

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        ImpedancePoint impedancePoint;
        impedancePoint.set_magnitude(5000);
        ImpedanceAtFrequency impedanceAtFrequency;
        impedanceAtFrequency.set_frequency(1e5);
        impedanceAtFrequency.set_impedance(impedancePoint);
        inputs.get_mutable_design_requirements().set_minimum_impedance(std::vector<ImpedanceAtFrequency>{impedanceAtFrequency});
        inputs.get_mutable_design_requirements().get_mutable_magnetizing_inductance().set_minimum(magnetizingInductance);
        inputs.get_mutable_design_requirements().get_mutable_magnetizing_inductance().set_nominal(std::nullopt);
        inputs.get_mutable_design_requirements().get_mutable_magnetizing_inductance().set_maximum(std::nullopt);

        OpenMagnetics::Mas masMagnetic;
        inputs.process();

        MagneticAdviser magneticAdviser;
        magneticAdviser.set_application(MAS::Application::INTERFERENCE_SUPPRESSION);
        magneticAdviser.set_core_mode(CoreAdviser::CoreAdviserModes::STANDARD_CORES);
        auto masMagnetics = magneticAdviser.get_advised_magnetic(inputs, 1);
        REQUIRE(masMagnetics.size() > 0);

        for (auto [masMagnetic, scoring] : masMagnetics) {
            auto impedance = Impedance().calculate_impedance(masMagnetic.get_mutable_magnetic(), 1e6);
            REQUIRE(abs(impedance) >= 5000);

            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            if (plot) {
                auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
                // {
                //     auto outFile = outputFilePath;
                //     std::string filename = "Test_MagneticAdviser_Filter_" + std::to_string(scoring) + ".svg";
                //     outFile.append(filename);
                //     settings.set_painter_include_fringing(false);
                //     Painter painter(outFile, true);

                //     painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
                //     painter.paint_core(masMagnetic.get_mutable_magnetic());
                //     painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
                //     painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
                //     painter.export_svg();
                // }
                {
                    auto impedanceSweep = Sweeper().sweep_impedance_over_frequency(masMagnetic.get_mutable_magnetic(), 10000, 400e6, 1000);

                    auto outFile = outputFilePath;

                    outFile.append("Test_MagneticAdviser_Filter_Sweep_Impedance_Over_Frequency_" + std::to_string(scoring) + ".svg");
                    std::filesystem::remove(outFile);
                    Painter painter(outFile, false, true);
                    #ifdef ENABLE_MATPLOTPP
                    painter.paint_curve(impedanceSweep);
                    #else
                        INFO("matplotplusplus disabled — skipping AdvancedPainter call");
                    #endif

                    #ifdef ENABLE_MATPLOTPP
                    painter.export_svg();
                    #else
                        INFO("matplotplusplus disabled — skipping AdvancedPainter call");
                    #endif

                }
            }
        }
        settings.reset();
    }

    TEST_CASE("Test_MagneticAdviser_Planar", "[adviser][magnetic-adviser]") {
        clear_databases();
        settings.reset();
        // settings.set_use_concentric_cores(false);
        // settings.set_use_toroidal_cores(true);
        // settings.set_use_only_cores_in_stock(false);
        settings.set_coil_maximum_layers_planar(17);

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {4, 1};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 507026;
        double magnetizingInductance = 10e-6;
        double temperature = 25;
        WaveformLabel waveShape = WaveformLabel::TRIANGULAR;
        double peakToPeak = 200;
        double dutyCycle = 0.5;
        double dcCurrent = 0;

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(frequency,
                                                                          magnetizingInductance,
                                                                          temperature,
                                                                          waveShape,
                                                                          peakToPeak,
                                                                          dutyCycle,
                                                                          dcCurrent,
                                                                          turnsRatios);
        inputs.get_mutable_design_requirements().get_mutable_magnetizing_inductance().set_minimum(magnetizingInductance);
        inputs.get_mutable_design_requirements().get_mutable_magnetizing_inductance().set_nominal(std::nullopt);
        inputs.get_mutable_design_requirements().get_mutable_magnetizing_inductance().set_maximum(std::nullopt);

        inputs.get_mutable_design_requirements().set_wiring_technology(MAS::WiringTechnology::PRINTED);

        OpenMagnetics::Mas masMagnetic;
        inputs.process();

        MagneticAdviser magneticAdviser;
        magneticAdviser.set_core_mode(CoreAdviser::CoreAdviserModes::STANDARD_CORES);
        auto masMagnetics = magneticAdviser.get_advised_magnetic(inputs, 10);
        REQUIRE(masMagnetics.size() > 0);

        for (auto [masMagnetic, scoring] : masMagnetics) {

            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            if (plot) {
                auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
                {
                    auto outFile = outputFilePath;
                    std::string filename = "Test_MagneticAdviser_Filter_Planar_" + std::to_string(scoring) + ".svg";
                    outFile.append(filename);
                    settings.set_painter_include_fringing(false);
                    Painter painter(outFile, false, false);

                    // painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
                    painter.paint_core(masMagnetic.get_mutable_magnetic());
                    painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
                    painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
                    #ifdef ENABLE_MATPLOTPP
                    painter.export_svg();
                    #else
                        INFO("matplotplusplus disabled — skipping AdvancedPainter call");
                    #endif

                }
            }
        }
        settings.reset();
    }

    TEST_CASE("Test_MagneticAdviser_High_Current", "[adviser][magnetic-adviser]") {

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {24, 78, 76};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 507026;
        double magnetizingInductance = 100e-6;
        double temperature = 25;
        WaveformLabel waveShape = WaveformLabel::TRIANGULAR;
        double peakToPeak = 100;
        double dutyCycle = 0.5;
        double dcCurrent = 0;

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        OpenMagnetics::Mas masMagnetic;
        inputs.process();

        MagneticAdviser magneticAdviser;
        auto masMagnetics = magneticAdviser.get_advised_magnetic(inputs, 5);

        REQUIRE(masMagnetics.size() > 0);
        for (auto [masMagnetic, scoring] : masMagnetics) {
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            if (plot) {
                auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
                auto outFile = outputFilePath;
                std::string filename = "Test_MagneticAdviser_High_Current_" + std::to_string(scoring) + ".svg";
                outFile.append(filename);
                settings.set_painter_include_fringing(false);
                Painter painter(outFile, true);

                painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
                painter.paint_core(masMagnetic.get_mutable_magnetic());
                painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
                painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
                #ifdef ENABLE_MATPLOTPP
                painter.export_svg();
                #else
                    INFO("matplotplusplus disabled — skipping AdvancedPainter call");
                #endif

            }
        }
        settings.reset();
    }

    TEST_CASE("Test_MagneticAdviser", "[adviser][magnetic-adviser]") {

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {24, 78, 76};
        std::vector<int64_t> numberParallels = {1, 1, 1};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 507026;
        double magnetizingInductance = 100e-6;
        double temperature = 25;
        WaveformLabel waveShape = WaveformLabel::TRIANGULAR;
        double peakToPeak = 1;
        double dutyCycle = 0.5;
        double dcCurrent = 0;

        settings.set_coil_allow_margin_tape(true);
        settings.set_coil_allow_insulated_wire(false);
        settings.set_coil_fill_sections_with_margin_tape(true);

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);


        auto requirements = inputs.get_mutable_design_requirements();
        auto insulationRequirements = requirements.get_insulation().value();
        auto standards = std::vector<InsulationStandards>{InsulationStandards::IEC_606641};
        insulationRequirements.set_standards(standards);
        requirements.set_insulation(insulationRequirements);
        inputs.set_design_requirements(requirements);

        OpenMagnetics::Mas masMagnetic;
        inputs.process();

        MagneticAdviser magneticAdviser;
        auto masMagnetics = magneticAdviser.get_advised_magnetic(inputs, 4);
        REQUIRE(masMagnetics.size() > 0);

        auto scorings = magneticAdviser.get_scorings();
        for (auto [name, values] : scorings) {
            double scoringTotal = 0;
            for (auto [key, value] : values) {
                scoringTotal += value;
            }
        }

        for (auto masMagneticWithScoring : masMagnetics) {
            auto masMagnetic = masMagneticWithScoring.first;
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            if (plot) {
                auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
                auto outFile = outputFilePath;
                std::string filename = "MagneticAdviser" + masMagnetic.get_magnetic().get_manufacturer_info().value().get_reference().value() + ".svg";
                outFile.append(filename);
                Painter painter(outFile, true);

                painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
                painter.paint_core(masMagnetic.get_mutable_magnetic());
                painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
                painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
                #ifdef ENABLE_MATPLOTPP
                painter.export_svg();
                #else
                    INFO("matplotplusplus disabled — skipping AdvancedPainter call");
                #endif

            }
        }
    }

    TEST_CASE("Test_MagneticAdviser_No_Insulation_Requirements", "[adviser][magnetic-adviser]") {

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {24, 78, 76};
        std::vector<int64_t> numberParallels = {1, 1, 1};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 507026;
        double magnetizingInductance = 100e-6;
        double temperature = 25;
        WaveformLabel waveShape = WaveformLabel::TRIANGULAR;
        double peakToPeak = 1;
        double dutyCycle = 0.5;
        double dcCurrent = 0;

        settings.set_coil_allow_margin_tape(true);
        settings.set_coil_allow_insulated_wire(false);
        settings.set_coil_fill_sections_with_margin_tape(true);

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);


        auto requirements = inputs.get_mutable_design_requirements();
        requirements.set_insulation(std::nullopt);
        inputs.set_design_requirements(requirements);

        OpenMagnetics::Mas masMagnetic;
        inputs.process();

        MagneticAdviser magneticAdviser;
        auto masMagnetics = magneticAdviser.get_advised_magnetic(inputs, 1);
        REQUIRE(masMagnetics.size() > 0);

        for (auto masMagneticWithScoring : masMagnetics) {
            auto masMagnetic = masMagneticWithScoring.first;
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            if (plot) {
                auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
                auto outFile = outputFilePath;
                std::string filename = "Test_MagneticAdviser_No_Insulation_Requirements_" + std::to_string(OpenMagnetics::TestUtils::randomInt(0, RAND_MAX)) + ".svg";
                outFile.append(filename);
                Painter painter(outFile, true);

                painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
                painter.paint_core(masMagnetic.get_mutable_magnetic());
                painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
                painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
                #ifdef ENABLE_MATPLOTPP
                painter.export_svg();
                #else
                    INFO("matplotplusplus disabled — skipping AdvancedPainter call");
                #endif

            }
        }
    }

    TEST_CASE("MagneticAdviserJsonHV", "[adviser][magnetic-adviser][smoke-test]") {
        clear_databases();
        settings.reset();
        settings.set_use_only_cores_in_stock(true);
        settings.set_use_toroidal_cores(true);
        settings.set_use_concentric_cores(true);
        settings.set_coil_fill_sections_with_margin_tape(true);
        std::string masString = R"({"inputs": {"designRequirements": {"insulation": {"altitude": {"maximum": 1000 }, "cti": "Group IIIB", "pollutionDegree": "P2", "overvoltageCategory": "OVC-II", "insulationType": "Basic", "mainSupplyVoltage": {"nominal": null, "minimum": 100, "maximum": 815 }, "standards": ["IEC 60664-1"] }, "leakageInductance": [{"maximum": 0.0000027 }, {"maximum": 7.2e-7 } ], "magnetizingInductance": {"maximum": null, "minimum": null, "nominal": 0.0007 }, "market": "Commercial", "maximumDimensions": {"width": null, "height": 0.02, "depth": null }, "maximumWeight": null, "name": "My Design Requirements", "operatingTemperature": null, "strayCapacitance": [{"maximum": 5e-11 }, {"maximum": 5e-11 } ], "terminalType": ["Pin", "Pin", "Pin"], "topology": "Flyback Converter", "turnsRatios": [{"nominal": 3.6 }, {"nominal": 7 } ] }, "operatingPoints": [{"name": "Operating Point No. 1", "conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"name": "Primary winding excitation", "frequency": 42000, "current": {"waveform": {"ancillaryLabel": null, "data": [0, 0, 1.1, 0, 0 ], "numberPeriods": null, "time": [0, 0, 0.000007, 0.000007, 0.00002380952380952381 ] }, "processed": {"dutyCycle": 0.294, "peakToPeak": 1.1, "offset": 0, "label": "Flyback Primary", "acEffectiveFrequency": 320324.70197566116, "effectiveFrequency": 299906.6380380921, "peak": 1.0815263605442182, "rms": 0.34251379452955516, "thd": 1.070755363168311 }, "harmonics": {"amplitudes": [0.16053906914328236, 0.29202196202762337, 0.21752279286410783, 0.1310845319633303, 0.07713860380309742, 0.06991594662329052, 0.06385479128189217, 0.04886550825820873, 0.041656942163512685, 0.04106933856014564, 0.036071164493133245, 0.03055299653431604, 0.02975135660550652, 0.028391215457215468, 0.02482068925428346, 0.023301839916889198, 0.023120564125526488, 0.021217296242203084, 0.019422484452150467, 0.019302597014669035, 0.018564086064230358, 0.016990298638197066, 0.016530128240448857, 0.016408746965119517, 0.015342013775419296, 0.014571250724490458, 0.014600261585603161, 0.014080719038131833, 0.013228104677648548, 0.013112646237623088, 0.012998033401152062, 0.012300975273665412, 0.011953654437324483, 0.012016309039770377, 0.011608286047010533, 0.011111174591424324, 0.011137874398755002, 0.01101797292674616, 0.010530630511366074, 0.01040049906973778, 0.01046208177619542, 0.010124703775041825, 0.009838808475080322, 0.009931053977331523, 0.0098034018119726, 0.009456464631666566, 0.009454698580407995, 0.009502890216301076, 0.009218272478811467, 0.00907613407579727, 0.009199516076645192, 0.00906471876755199, 0.008824889331959554, 0.008907818118804589, 0.008937898630441313, 0.00869876368067147, 0.008666288051213203, 0.008803433740422266, 0.008663604885115571, 0.008515486625290338, 0.008660125560530002, 0.008670306514975606, 0.008475927879718145, 0.008536460912837853 ], "frequencies": [0, 42000, 84000, 126000, 168000, 210000, 252000, 294000, 336000, 378000, 420000, 462000, 504000, 546000, 588000, 630000, 672000, 714000, 756000, 798000, 840000, 882000, 924000, 966000, 1008000, 1050000, 1092000, 1134000, 1176000, 1218000, 1260000, 1302000, 1344000, 1386000, 1428000, 1470000, 1512000, 1554000, 1596000, 1638000, 1680000, 1722000, 1764000, 1806000, 1848000, 1890000, 1932000, 1974000, 2016000, 2058000, 2100000, 2142000, 2184000, 2226000, 2268000, 2310000, 2352000, 2394000, 2436000, 2478000, 2520000, 2562000, 2604000, 2646000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-29.4, 70.6, 70.6, -29.4, -29.4 ], "numberPeriods": null, "time": [0, 0, 0.000007, 0.000007, 0.00002380952380952381 ] }, "processed": {"dutyCycle": 0.294, "peakToPeak": 100, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 273993.1760953399, "effectiveFrequency": 273985.04815560044, "peak": 70.6, "rms": 45.33538904652732, "thd": 0.7943280084778311 }, "harmonics": {"amplitudes": [0.4937500000000009, 50.19273146886441, 30.88945734421223, 8.60726309175931, 7.514577540357868, 12.576487416046705, 7.890222788854489, 0.6723406687769768, 6.659332230305735, 6.771471403422372, 2.1663923211140617, 3.134209960512532, 5.356728264586365, 3.4347141143172375, 0.6805380572551197, 3.777520625726094, 3.7722086912079633, 1.0283716619016732, 2.1769829843489013, 3.474177141872247, 2.1027712520453545, 0.6945565902694819, 2.7474727308098412, 2.608727438661375, 0.5486770349329464, 1.7724349572072589, 2.6198081625264167, 1.4623766283543616, 0.71496628279737, 2.2318766983434566, 1.995656440067472, 0.2773827898281461, 1.562500000000001, 2.1411881515502134, 1.0841275629577949, 0.7426431917434089, 1.9342807953981735, 1.620450382955772, 0.09545259242683589, 1.4471144811333827, 1.8430944173488684, 0.8315259963159228, 0.7788650462348548, 1.7521633971823452, 1.3695427043271637, 0.042929496561759636, 1.388303796831737, 1.6472150875774871, 0.6472086912079588, 0.8254655623504514, 1.6415464278026473, 1.1917877639027972, 0.16004318769545572, 1.3696996040394052, 1.5166136111631239, 0.5023333497836172, 0.8850851521101957, 1.581550490460967, 1.0607923422947243, 0.2691521002455623, 1.3846695378376832, 1.4323337135341552, 0.38011439553529425, 0.9615889729028446 ], "frequencies": [0, 42000, 84000, 126000, 168000, 210000, 252000, 294000, 336000, 378000, 420000, 462000, 504000, 546000, 588000, 630000, 672000, 714000, 756000, 798000, 840000, 882000, 924000, 966000, 1008000, 1050000, 1092000, 1134000, 1176000, 1218000, 1260000, 1302000, 1344000, 1386000, 1428000, 1470000, 1512000, 1554000, 1596000, 1638000, 1680000, 1722000, 1764000, 1806000, 1848000, 1890000, 1932000, 1974000, 2016000, 2058000, 2100000, 2142000, 2184000, 2226000, 2268000, 2310000, 2352000, 2394000, 2436000, 2478000, 2520000, 2562000, 2604000, 2646000 ] } } }, {"name": "Primary winding excitation", "frequency": 42000, "current": {"waveform": {"ancillaryLabel": null, "data": [0, 0, 3.96, 0, 0 ], "numberPeriods": null, "time": [0, 0.000011904761904761905, 0.000011904761904761905, 0.00002380952380952381, 0.00002380952380952381 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 3.96, "offset": 0, "label": "Flyback Secondary", "acEffectiveFrequency": 273605.22552591236, "effectiveFrequency": 239619.35667760254, "peak": 3.9599999999999986, "rms": 1.635596311125932, "thd": 0.6765065136547783 }, "harmonics": {"amplitudes": [1.0054687500000032, 1.5109819803933975, 0.6305067526445318, 0.43631079159665215, 0.3156335707813906, 0.25867345007829645, 0.2108457708927511, 0.1845098972167098, 0.15858039332900484, 0.14374429235357594, 0.12732519999086678, 0.11798991339580255, 0.1065764142034965, 0.10027664365108102, 0.09183268563443873, 0.08737513814977015, 0.08084358345172529, 0.0775829517661027, 0.07235909588115942, 0.06991721285428486, 0.06562943182065245, 0.06377074179801741, 0.060177659717708575, 0.05874824283059186, 0.05568602880656547, 0.054581360149729106, 0.051934759310085545, 0.0510815825974439, 0.0487670673260459, 0.048112679051860716, 0.04606818176826717, 0.04557382535330518, 0.043752232085917345, 0.04338887874386873, 0.04175376612752825, 0.04149933709574444, 0.04002209784539918, 0.03985958700992044, 0.038517442612831645, 0.03843361846308103, 0.03720821487910577, 0.037192705080064635, 0.03606910115775439, 0.036113735798427, 0.03507966215426201, 0.03517799562893598, 0.03422330356713118, 0.0343702622545533, 0.03348650869654576, 0.03367812888504953, 0.03285826032368028, 0.03309149205770409, 0.032329601799393176, 0.032602161771072315, 0.0318933023014623, 0.032203563952980456, 0.031543601457069734, 0.03189051394290526, 0.03127601564366044, 0.03165904576198601, 0.03108719333289207, 0.0315062863179612, 0.03097481051580278, 0.03143036691725387 ], "frequencies": [0, 42000, 84000, 126000, 168000, 210000, 252000, 294000, 336000, 378000, 420000, 462000, 504000, 546000, 588000, 630000, 672000, 714000, 756000, 798000, 840000, 882000, 924000, 966000, 1008000, 1050000, 1092000, 1134000, 1176000, 1218000, 1260000, 1302000, 1344000, 1386000, 1428000, 1470000, 1512000, 1554000, 1596000, 1638000, 1680000, 1722000, 1764000, 1806000, 1848000, 1890000, 1932000, 1974000, 2016000, 2058000, 2100000, 2142000, 2184000, 2226000, 2268000, 2310000, 2352000, 2394000, 2436000, 2478000, 2520000, 2562000, 2604000, 2646000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-13.89, 13.89, 13.89, -13.89, -13.89 ], "numberPeriods": null, "time": [0, 0, 0.000011904761904761905, 0.000011904761904761905, 0.00002380952380952381 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 27.78, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 248423.92512497256, "effectiveFrequency": 248408.7565140517, "peak": 13.89, "rms": 13.890000000000022, "thd": 0.48331514845248535 }, "harmonics": {"amplitudes": [0.21703125, 17.681745968226167, 0.4340625, 5.884441743008607, 0.4340625, 3.5192857754085134, 0.4340625, 2.50156382659689, 0.4340625, 1.932968090534883, 0.4340625, 1.56850033166751, 0.4340625, 1.313925940874188, 0.4340625, 1.1252647178556892, 0.4340625, 0.9792293094780014, 0.4340625, 0.8623340820515439, 0.4340625, 0.7662274695502621, 0.4340625, 0.6854595927802318, 0.4340625, 0.6163213952979067, 0.4340625, 0.5561996920846202, 0.4340625, 0.5031990666519233, 0.4340625, 0.45591010107099894, 0.4340625, 0.41326185461487136, 0.4340625, 0.3744248874701925, 0.4340625, 0.33874569976854113, 0.4340625, 0.30570130348173624, 0.4340625, 0.2748670467095772, 0.4340625, 0.24589336899764186, 0.4340625, 0.21848870156913297, 0.4340625, 0.1924066733732579, 0.4340625, 0.16743638267205946, 0.4340625, 0.1433948809785266, 0.4340625, 0.12012127132032469, 0.4340625, 0.09747199388796957, 0.4340625, 0.07531698847858959, 0.4340625, 0.05353650312310587, 0.4340625, 0.03201837355771753, 0.4340625, 0.0106556362841701 ], "frequencies": [0, 42000, 84000, 126000, 168000, 210000, 252000, 294000, 336000, 378000, 420000, 462000, 504000, 546000, 588000, 630000, 672000, 714000, 756000, 798000, 840000, 882000, 924000, 966000, 1008000, 1050000, 1092000, 1134000, 1176000, 1218000, 1260000, 1302000, 1344000, 1386000, 1428000, 1470000, 1512000, 1554000, 1596000, 1638000, 1680000, 1722000, 1764000, 1806000, 1848000, 1890000, 1932000, 1974000, 2016000, 2058000, 2100000, 2142000, 2184000, 2226000, 2268000, 2310000, 2352000, 2394000, 2436000, 2478000, 2520000, 2562000, 2604000, 2646000 ] } } }, {"name": "Primary winding excitation", "frequency": 42000, "current": {"waveform": {"ancillaryLabel": null, "data": [0, 0, 0.08, 0, 0 ], "numberPeriods": null, "time": [0, 0.000011904761904761905, 0.000011904761904761905, 0.00002380952380952381, 0.00002380952380952381 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 0.08, "offset": 0, "label": "Flyback Secondary", "acEffectiveFrequency": 273605.2255259123, "effectiveFrequency": 239619.35667760245, "peak": 0.07999999999999997, "rms": 0.03304234971971578, "thd": 0.6765065136547783 }, "harmonics": {"amplitudes": [0.020312500000000067, 0.030524888492795912, 0.012737510154434986, 0.008814359426194993, 0.006376435773361427, 0.00522572626420801, 0.004259510523085881, 0.0037274726710446416, 0.0032036443096768643, 0.0029039250980520386, 0.0025722262624417523, 0.002383634614056616, 0.002153058872797908, 0.0020257907808299194, 0.0018552057703927027, 0.0017651543060559615, 0.0016332037060954627, 0.0015673323589111663, 0.0014617999167910994, 0.0014124689465512088, 0.001325847107487929, 0.0012882978141013626, 0.0012157102973274464, 0.001186833188496806, 0.001124970278920515, 0.0011026537403985668, 0.0010491870567694052, 0.0010319511635847236, 0.0009851932793140594, 0.000971973314179004, 0.0009306703387528712, 0.0009206833404708095, 0.0008838834764831795, 0.0008765430049266404, 0.0008435104268187535, 0.0008383704463786768, 0.000808527229199984, 0.0008052441820185946, 0.0007781301537945788, 0.0007764367366278998, 0.0007516811086688038, 0.0007513677793952459, 0.0007286687102576646, 0.0007295704201702419, 0.0007086800435204446, 0.0007106665783623426, 0.000691379870043054, 0.0006943487324152177, 0.0006764951251827419, 0.00068036624010201, 0.0006638032388622281, 0.0006685149910647295, 0.000653123268674609, 0.0006586295307287333, 0.0006443091374032789, 0.0006505770495551611, 0.0006372444738801973, 0.0006442528069273796, 0.0006318386998719282, 0.0006395766820603246, 0.0006280241077351936, 0.0006364906326860852, 0.0006257537477939952, 0.0006349569074192705 ], "frequencies": [0, 42000, 84000, 126000, 168000, 210000, 252000, 294000, 336000, 378000, 420000, 462000, 504000, 546000, 588000, 630000, 672000, 714000, 756000, 798000, 840000, 882000, 924000, 966000, 1008000, 1050000, 1092000, 1134000, 1176000, 1218000, 1260000, 1302000, 1344000, 1386000, 1428000, 1470000, 1512000, 1554000, 1596000, 1638000, 1680000, 1722000, 1764000, 1806000, 1848000, 1890000, 1932000, 1974000, 2016000, 2058000, 2100000, 2142000, 2184000, 2226000, 2268000, 2310000, 2352000, 2394000, 2436000, 2478000, 2520000, 2562000, 2604000, 2646000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-7.14, 7.14, 7.14, -7.14, -7.14 ], "numberPeriods": null, "time": [0, 0, 0.000011904761904761905, 0.000011904761904761905, 0.00002380952380952381 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 14.28, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 248423.9251249721, "effectiveFrequency": 248408.75651405123, "peak": 7.14, "rms": 7.139999999999993, "thd": 0.4833151484524845 }, "harmonics": {"amplitudes": [0.1115625, 9.08910483895859, 0.223125, 3.0248318246998886, 0.223125, 1.8090497074454124, 0.223125, 1.2859010598921377, 0.223125, 0.9936207463224668, 0.223125, 0.8062701488917218, 0.223125, 0.6754090149634054, 0.223125, 0.5784298117703113, 0.223125, 0.5033619344616945, 0.223125, 0.44327324304161464, 0.223125, 0.3938707078897675, 0.223125, 0.3523528792261233, 0.223125, 0.3168131578421202, 0.223125, 0.285908265045658, 0.223125, 0.2586638830737746, 0.223125, 0.23435551631727392, 0.223125, 0.21243265960764457, 0.223125, 0.19246894863478564, 0.223125, 0.1741284590602866, 0.223125, 0.1571423547055144, 0.223125, 0.14129234798462062, 0.223125, 0.1263987512342087, 0.223125, 0.1123116867677183, 0.223125, 0.09890451028690123, 0.223125, 0.08606881009924447, 0.223125, 0.07371054356995504, 0.223125, 0.0617470034000806, 0.223125, 0.05010439426638574, 0.223125, 0.0387158601682599, 0.223125, 0.027519843938011213, 0.223125, 0.016458688783449027, 0.223125, 0.005477411308063118 ], "frequencies": [0, 42000, 84000, 126000, 168000, 210000, 252000, 294000, 336000, 378000, 420000, 462000, 504000, 546000, 588000, 630000, 672000, 714000, 756000, 798000, 840000, 882000, 924000, 966000, 1008000, 1050000, 1092000, 1134000, 1176000, 1218000, 1260000, 1302000, 1344000, 1386000, 1428000, 1470000, 1512000, 1554000, 1596000, 1638000, 1680000, 1722000, 1764000, 1806000, 1848000, 1890000, 1932000, 1974000, 2016000, 2058000, 2100000, 2142000, 2184000, 2226000, 2268000, 2310000, 2352000, 2394000, 2436000, 2478000, 2520000, 2562000, 2604000, 2646000 ] } } } ] } ] }, "magnetic": {"coil": {"functionalDescription": [{"name": "Primary"}, {"name": "Secondary"}, {"name": "Tertiary"} ] } } })";
        json masJson = json::parse(masString);

        // if (masJson["inputs"]["designRequirements"]["insulation"]["cti"] == "GROUP_IIIa") {
        //     masJson["inputs"]["designRequirements"]["insulation"]["cti"] = "GROUP_IIIA";
        // }
        // if (masJson["inputs"]["designRequirements"]["insulation"]["cti"] == "GROUP_IIIb") {
        //     masJson["inputs"]["designRequirements"]["insulation"]["cti"] = "GROUP_IIIB";
        // }

        OpenMagnetics::Inputs inputs(masJson["inputs"]);
        OpenMagnetics::Mas masMagnetic;

        MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 2);
        REQUIRE(masMagnetics.size() > 0);

        for (auto [masMagnetic, scoring] : masMagnetics) {
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "MagneticAdviser_MagneticAdviserJsonHV_" + std::to_string(scoring) + ".svg";
            outFile.append(filename);
            Painter painter(outFile, true);

            painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
            painter.paint_core(masMagnetic.get_mutable_magnetic());
            painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            #ifdef ENABLE_MATPLOTPP
            painter.export_svg();
            #else
                INFO("matplotplusplus disabled — skipping AdvancedPainter call");
            #endif

        }
    }

    TEST_CASE("MagneticAdviserJsonLV", "[adviser][magnetic-adviser][smoke-test]") {
        clear_databases();
        settings.reset();
        settings.set_use_only_cores_in_stock(true);
        settings.set_coil_fill_sections_with_margin_tape(true);
        std::string masString = R"({"inputs": {"designRequirements": {"insulation": {"altitude": {"maximum": 1000 }, "cti": "Group IIIB", "pollutionDegree": "P2", "overvoltageCategory": "OVC-II", "insulationType": "Basic", "mainSupplyVoltage": {"nominal": 120, "minimum": 20, "maximum": 450 }, "standards": ["IEC 60664-1"] }, "leakageInductance": [{"maximum": 0.00000135 }, {"maximum": 0.00000135 } ], "magnetizingInductance": {"maximum": 0.000088, "minimum": 0.000036, "nominal": 0.000039999999999999996 }, "market": "Commercial", "maximumDimensions": {"width": null, "height": 0.02, "depth": null }, "maximumWeight": null, "name": "My Design Requirements", "operatingTemperature": null, "strayCapacitance": [{"maximum": 5e-11 }, {"maximum": 5e-11 } ], "terminalType": ["Pin", "Pin", "Pin"], "topology": "Flyback Converter", "turnsRatios": [{"nominal": 1.2 }, {"nominal": 1.2 } ] }, "operatingPoints": [{"name": "Operating Point No. 1", "conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"name": "Primary winding excitation", "frequency": 36000, "current": {"waveform": {"ancillaryLabel": null, "data": [0, 0, 5.5, 0, 0 ], "numberPeriods": null, "time": [0, 0, 0.000011111111111111112, 0.000011111111111111112, 0.00002777777777777778 ] }, "processed": {"dutyCycle": 0.4, "peakToPeak": 5.5, "offset": 0, "label": "Flyback Primary", "acEffectiveFrequency": 247671.20334838802, "effectiveFrequency": 224591.09893738205, "peak": 5.478515625000006, "rms": 2.025897603604934, "thd": 0.8123116860415294 }, "harmonics": {"amplitudes": [1.1128234863281259, 1.8581203560604207, 1.051281221517565, 0.5159082238421047, 0.4749772671818343, 0.35101798151542746, 0.28680520988943464, 0.26606329358167974, 0.2109058168232097, 0.2038324242201338, 0.17683391011628416, 0.15977271740890614, 0.15330090420784062, 0.13359870506402888, 0.1311678604251054, 0.11938227872282572, 0.11196939263234992, 0.10881687126769939, 0.09881440376099933, 0.09781517962337143, 0.09113803697574208, 0.08720137281111257, 0.08534603089826898, 0.07933100843010468, 0.07894777491614259, 0.0746052794656816, 0.07229734773884018, 0.07109041486606842, 0.0670974294019412, 0.06703716258861198, 0.06396085871081211, 0.0625474905371643, 0.06171647321147238, 0.05889230550248535, 0.059028047581821314, 0.05671822612277054, 0.05585484049699437, 0.055265364667603324, 0.0531817689812044, 0.053452567468547284, 0.05164394009764037, 0.051150939884780554, 0.05072953903613357, 0.049149522805038266, 0.049525208170729855, 0.04806384343218855, 0.04784167797419847, 0.047545488047249944, 0.04632889186784806, 0.046795775944523485, 0.045586400651502844, 0.045578916197115495, 0.045382141359988915, 0.04444264358178238, 0.04499725549277023, 0.04397809673849301, 0.04415534166549457, 0.044042673430388084, 0.04332652563384539, 0.04397294643650677, 0.043102730085842, 0.043452676225114936, 0.04341597260904029, 0.04289100286968138 ], "frequencies": [0, 36000, 72000, 108000, 144000, 180000, 216000, 252000, 288000, 324000, 360000, 396000, 432000, 468000, 504000, 540000, 576000, 612000, 648000, 684000, 720000, 756000, 792000, 828000, 864000, 900000, 936000, 972000, 1008000, 1044000, 1080000, 1116000, 1152000, 1188000, 1224000, 1260000, 1296000, 1332000, 1368000, 1404000, 1440000, 1476000, 1512000, 1548000, 1584000, 1620000, 1656000, 1692000, 1728000, 1764000, 1800000, 1836000, 1872000, 1908000, 1944000, 1980000, 2016000, 2052000, 2088000, 2124000, 2160000, 2196000, 2232000, 2268000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-8, 12, 12, -8, -8 ], "numberPeriods": null, "time": [0, 0, 0.000011111111111111112, 0.000011111111111111112, 0.00002777777777777778 ] }, "processed": {"dutyCycle": 0.4, "peakToPeak": 20, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 217445.61780293446, "effectiveFrequency": 217445.06394458748, "peak": 12, "rms": 9.79157801378307, "thd": 0.5579294461931013 }, "harmonics": {"amplitudes": [0.03125, 12.090982167347164, 3.7938629699811433, 2.4460154592542187, 3.050934294867692, 0.06265085868474476, 2.0052584267980813, 1.1245773748206938, 0.8899247078195547, 1.3746157728063957, 0.06310658027096225, 1.0931619934351904, 0.682943044655975, 0.5329821111108219, 0.8998037891938968, 0.06387675236155328, 0.7544417382415929, 0.5036912732460003, 0.37575767434186663, 0.6781566555968945, 0.06497787099509279, 0.579613707579126, 0.4082105683349906, 0.28790910341678466, 0.5516771452281819, 0.06643416431010025, 0.4742269329016237, 0.35024440150573616, 0.23220844850577863, 0.4713901788294872, 0.0682788501833682, 0.40480954492681936, 0.3125, 0.19399945657286152, 0.41719035005974064, 0.07055595095604376, 0.35652867811806305, 0.28709106236926035, 0.16634675530328177, 0.3793488450042862, 0.07332285477181877, 0.32183350571658953, 0.2699539504381825, 0.14554626367323475, 0.3526337777044684, 0.0766539126369276, 0.29650814745710025, 0.25883312123614344, 0.12944173824159266, 0.33403609127672834, 0.08064551519955634, 0.27805045130481554, 0.25243555700995646, 0.1166943265232067, 0.3217665194149353, 0.08542334196565987, 0.2649247078195551, 0.250039428516803, 0.10643002348369442, 0.31477312120987916, 0.09115288771096393, 0.25618941041585475, 0.2513050615401854, 0.098055076230366 ], "frequencies": [0, 36000, 72000, 108000, 144000, 180000, 216000, 252000, 288000, 324000, 360000, 396000, 432000, 468000, 504000, 540000, 576000, 612000, 648000, 684000, 720000, 756000, 792000, 828000, 864000, 900000, 936000, 972000, 1008000, 1044000, 1080000, 1116000, 1152000, 1188000, 1224000, 1260000, 1296000, 1332000, 1368000, 1404000, 1440000, 1476000, 1512000, 1548000, 1584000, 1620000, 1656000, 1692000, 1728000, 1764000, 1800000, 1836000, 1872000, 1908000, 1944000, 1980000, 2016000, 2052000, 2088000, 2124000, 2160000, 2196000, 2232000, 2268000 ] } } }, {"name": "Primary winding excitation", "frequency": 36000, "current": {"waveform": {"ancillaryLabel": null, "data": [0, 0, 6.6, 0, 0 ], "numberPeriods": null, "time": [0, 0.00001388888888888889, 0.00001388888888888889, 0.00002777777777777778, 0.00002777777777777778 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 6.6, "offset": 0, "label": "Flyback Secondary", "acEffectiveFrequency": 234518.76473649617, "effectiveFrequency": 205388.02000937346, "peak": 6.599999999999992, "rms": 2.7259938518765505, "thd": 0.6765065136547781 }, "harmonics": {"amplitudes": [1.6757812500000036, 2.5183033006556617, 1.0508445877408874, 0.7271846526610849, 0.5260559513023173, 0.43112241679715924, 0.35140961815458516, 0.30751649536118214, 0.2643006555483416, 0.23957382058929272, 0.21220866665144478, 0.19664985565967041, 0.17762735700582763, 0.16712773941846806, 0.15305447605739808, 0.14562523024961643, 0.13473930575287565, 0.12930491961017093, 0.12059849313526576, 0.11652868809047452, 0.10938238636775419, 0.10628456966336206, 0.10029609952951436, 0.0979137380509861, 0.09281004801094254, 0.09096893358288156, 0.08655793218347599, 0.08513597099573951, 0.08127844554340995, 0.08018779841976777, 0.07678030294711186, 0.07595637558884157, 0.07292038680986246, 0.07231479790644763, 0.06958961021254718, 0.06916556182624073, 0.06670349640899874, 0.06643264501653386, 0.06419573768805274, 0.06405603077180153, 0.062013691465176275, 0.0619878418001075, 0.060115168596257366, 0.060189559664044844, 0.058466103590436726, 0.05862999271489314, 0.05703883927855185, 0.05728377042425536, 0.05581084782757627, 0.05613021480841572, 0.054763767206133816, 0.05515248676284005, 0.05388266966565528, 0.05433693628512042, 0.05315550383577058, 0.05367260658830061, 0.05257266909511635, 0.05315085657150875, 0.052126692739434176, 0.05276507626997662, 0.05181198888815346, 0.05251047719660189, 0.051624684193004804, 0.05238394486208969 ], "frequencies": [0, 36000, 72000, 108000, 144000, 180000, 216000, 252000, 288000, 324000, 360000, 396000, 432000, 468000, 504000, 540000, 576000, 612000, 648000, 684000, 720000, 756000, 792000, 828000, 864000, 900000, 936000, 972000, 1008000, 1044000, 1080000, 1116000, 1152000, 1188000, 1224000, 1260000, 1296000, 1332000, 1368000, 1404000, 1440000, 1476000, 1512000, 1548000, 1584000, 1620000, 1656000, 1692000, 1728000, 1764000, 1800000, 1836000, 1872000, 1908000, 1944000, 1980000, 2016000, 2052000, 2088000, 2124000, 2160000, 2196000, 2232000, 2268000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-8.335, 8.335, 8.335, -8.335, -8.335 ], "numberPeriods": null, "time": [0, 0, 0.00001388888888888889, 0.00001388888888888889, 0.00002777777777777778 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 16.67, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 212934.79296426187, "effectiveFrequency": 212921.79129775826, "peak": 8.335, "rms": 8.33500000000001, "thd": 0.48331514845248513 }, "harmonics": {"amplitudes": [0.130234375, 10.610320564806702, 0.26046875, 3.531088691718988, 0.26046875, 2.1118248335514735, 0.26046875, 1.5011183941457937, 0.26046875, 1.1599200168904427, 0.26046875, 0.9412131219905466, 0.26046875, 0.7884501596246478, 0.26046875, 0.6752398432920929, 0.26046875, 0.587608084557174, 0.26046875, 0.5174625323181874, 0.26046875, 0.4597916456948479, 0.26046875, 0.4113251048108877, 0.26046875, 0.36983720876947807, 0.26046875, 0.333759858425148, 0.26046875, 0.3019556674257579, 0.26046875, 0.2735788835440449, 0.26046875, 0.24798686524225744, 0.26046875, 0.2246818889174983, 0.26046875, 0.20327180760048874, 0.26046875, 0.18344279082219334, 0.26046875, 0.16494001687000165, 0.26046875, 0.14755372430492025, 0.26046875, 0.1311089508695984, 0.26046875, 0.1154578561962637, 0.26046875, 0.10047388405843227, 0.26046875, 0.0860472521926578, 0.26046875, 0.0720814108318868, 0.26046875, 0.058490213754947895, 0.26046875, 0.04519561547653339, 0.26046875, 0.032125756193742694, 0.26046875, 0.019213329273115587, 0.26046875, 0.006394148914942832 ], "frequencies": [0, 36000, 72000, 108000, 144000, 180000, 216000, 252000, 288000, 324000, 360000, 396000, 432000, 468000, 504000, 540000, 576000, 612000, 648000, 684000, 720000, 756000, 792000, 828000, 864000, 900000, 936000, 972000, 1008000, 1044000, 1080000, 1116000, 1152000, 1188000, 1224000, 1260000, 1296000, 1332000, 1368000, 1404000, 1440000, 1476000, 1512000, 1548000, 1584000, 1620000, 1656000, 1692000, 1728000, 1764000, 1800000, 1836000, 1872000, 1908000, 1944000, 1980000, 2016000, 2052000, 2088000, 2124000, 2160000, 2196000, 2232000, 2268000 ] } } }, {"name": "Primary winding excitation", "frequency": 36000, "current": {"waveform": {"ancillaryLabel": null, "data": [0, 0, 0.08, 0, 0 ], "numberPeriods": null, "time": [0, 0.00001388888888888889, 0.00001388888888888889, 0.00002777777777777778, 0.00002777777777777778 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 0.08, "offset": 0, "label": "Flyback Secondary", "acEffectiveFrequency": 234518.7647364962, "effectiveFrequency": 205388.02000937346, "peak": 0.0799999999999999, "rms": 0.033042349719715765, "thd": 0.6765065136547783 }, "harmonics": {"amplitudes": [0.020312500000000046, 0.030524888492795898, 0.012737510154435002, 0.008814359426194969, 0.006376435773361423, 0.005225726264207993, 0.004259510523085882, 0.003727472671044633, 0.0032036443096768678, 0.0029039250980520334, 0.002572226262441755, 0.002383634614056611, 0.0021530588727979106, 0.0020257907808299163, 0.001855205770392704, 0.0017651543060559591, 0.0016332037060954637, 0.0015673323589111632, 0.0014617999167910998, 0.0014124689465512055, 0.0013258471074879294, 0.0012882978141013591, 0.0012157102973274472, 0.0011868331884968021, 0.0011249702789205155, 0.0011026537403985642, 0.0010491870567694056, 0.0010319511635847216, 0.0009851932793140603, 0.0009719733141790018, 0.0009306703387528721, 0.0009206833404708108, 0.0008838834764831795, 0.0008765430049266376, 0.0008435104268187531, 0.0008383704463786751, 0.0008085272291999843, 0.0008052441820185928, 0.0007781301537945793, 0.0007764367366278984, 0.0007516811086688035, 0.0007513677793952425, 0.0007286687102576648, 0.0007295704201702394, 0.0007086800435204447, 0.0007106665783623409, 0.0006913798700430537, 0.0006943487324152148, 0.0006764951251827443, 0.0006803662401020098, 0.0006638032388622283, 0.0006685149910647278, 0.0006531232686746098, 0.0006586295307287323, 0.0006443091374032788, 0.0006505770495551569, 0.0006372444738801973, 0.0006442528069273779, 0.0006318386998719287, 0.0006395766820603229, 0.0006280241077351938, 0.0006364906326860843, 0.0006257537477939962, 0.0006349569074192687 ], "frequencies": [0, 36000, 72000, 108000, 144000, 180000, 216000, 252000, 288000, 324000, 360000, 396000, 432000, 468000, 504000, 540000, 576000, 612000, 648000, 684000, 720000, 756000, 792000, 828000, 864000, 900000, 936000, 972000, 1008000, 1044000, 1080000, 1116000, 1152000, 1188000, 1224000, 1260000, 1296000, 1332000, 1368000, 1404000, 1440000, 1476000, 1512000, 1548000, 1584000, 1620000, 1656000, 1692000, 1728000, 1764000, 1800000, 1836000, 1872000, 1908000, 1944000, 1980000, 2016000, 2052000, 2088000, 2124000, 2160000, 2196000, 2232000, 2268000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-8.335, 8.335, 8.335, -8.335, -8.335 ], "numberPeriods": null, "time": [0, 0, 0.00001388888888888889, 0.00001388888888888889, 0.00002777777777777778 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 16.67, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 212934.79296426187, "effectiveFrequency": 212921.79129775826, "peak": 8.335, "rms": 8.33500000000001, "thd": 0.48331514845248513 }, "harmonics": {"amplitudes": [0.130234375, 10.610320564806702, 0.26046875, 3.531088691718988, 0.26046875, 2.1118248335514735, 0.26046875, 1.5011183941457937, 0.26046875, 1.1599200168904427, 0.26046875, 0.9412131219905466, 0.26046875, 0.7884501596246478, 0.26046875, 0.6752398432920929, 0.26046875, 0.587608084557174, 0.26046875, 0.5174625323181874, 0.26046875, 0.4597916456948479, 0.26046875, 0.4113251048108877, 0.26046875, 0.36983720876947807, 0.26046875, 0.333759858425148, 0.26046875, 0.3019556674257579, 0.26046875, 0.2735788835440449, 0.26046875, 0.24798686524225744, 0.26046875, 0.2246818889174983, 0.26046875, 0.20327180760048874, 0.26046875, 0.18344279082219334, 0.26046875, 0.16494001687000165, 0.26046875, 0.14755372430492025, 0.26046875, 0.1311089508695984, 0.26046875, 0.1154578561962637, 0.26046875, 0.10047388405843227, 0.26046875, 0.0860472521926578, 0.26046875, 0.0720814108318868, 0.26046875, 0.058490213754947895, 0.26046875, 0.04519561547653339, 0.26046875, 0.032125756193742694, 0.26046875, 0.019213329273115587, 0.26046875, 0.006394148914942832 ], "frequencies": [0, 36000, 72000, 108000, 144000, 180000, 216000, 252000, 288000, 324000, 360000, 396000, 432000, 468000, 504000, 540000, 576000, 612000, 648000, 684000, 720000, 756000, 792000, 828000, 864000, 900000, 936000, 972000, 1008000, 1044000, 1080000, 1116000, 1152000, 1188000, 1224000, 1260000, 1296000, 1332000, 1368000, 1404000, 1440000, 1476000, 1512000, 1548000, 1584000, 1620000, 1656000, 1692000, 1728000, 1764000, 1800000, 1836000, 1872000, 1908000, 1944000, 1980000, 2016000, 2052000, 2088000, 2124000, 2160000, 2196000, 2232000, 2268000 ] } } } ] } ] }, "magnetic": {"coil": {"functionalDescription": [{"name": "Primary"}, {"name": "Secondary"}, {"name": "Tertiary"} ] } } })";
        json masJson = json::parse(masString);

        if (masJson["inputs"]["designRequirements"]["insulation"]["cti"] == "GROUP_IIIa") {
            masJson["inputs"]["designRequirements"]["insulation"]["cti"] = "GROUP_IIIA";
        }
        if (masJson["inputs"]["designRequirements"]["insulation"]["cti"] == "GROUP_IIIb") {
            masJson["inputs"]["designRequirements"]["insulation"]["cti"] = "GROUP_IIIB";
        }

        OpenMagnetics::Inputs inputs(masJson["inputs"]);
        OpenMagnetics::Mas masMagnetic;
        // inputs.process();

        MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 1);

        REQUIRE(masMagnetics.size() > 0);

        for (auto [masMagnetic, scoring] : masMagnetics) {
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            if (plot) {
                auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
                auto outFile = outputFilePath;
                std::string filename = "MagneticAdviser_MagneticAdviserJsonLV_" + std::to_string(scoring) + ".svg";
                outFile.append(filename);
                Painter painter(outFile, true);

                painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
                painter.paint_core(masMagnetic.get_mutable_magnetic());
                painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
                painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
                #ifdef ENABLE_MATPLOTPP
                painter.export_svg();
                #else
                    INFO("matplotplusplus disabled — skipping AdvancedPainter call");
                #endif

            }
        }
    }

    TEST_CASE("Test_MagneticAdviser_Random", "[adviser][magnetic-adviser][random]") {
        settings.reset();

        int count = 10;
        while (count > 0) {
            std::vector<double> turnsRatios;

            std::vector<int64_t> numberTurns;
            std::vector<IsolationSide> isolationSides;
            // for (size_t windingIndex = 0; windingIndex < OpenMagnetics::TestUtils::randomSize(1, 2 + 1 - 1); ++windingIndex)
            for (size_t windingIndex = 0; windingIndex < OpenMagnetics::TestUtils::randomSize(1, 4 + 1 - 1); ++windingIndex)
            {
                numberTurns.push_back(OpenMagnetics::TestUtils::randomInt64(1, 100 + 1 - 1));
                isolationSides.push_back(get_isolation_side_from_index(static_cast<size_t>(OpenMagnetics::TestUtils::randomInt64(1, 10 + 1 - 1))));
            }
            for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
                turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
            }

            double frequency = OpenMagnetics::TestUtils::randomInt64(10000, 1000000 + 10000 - 1);
            double magnetizingInductance = double(OpenMagnetics::TestUtils::randomInt(0, 10000 - 1)) *  1e-6;
            double temperature = 25;
            double peakToPeak = OpenMagnetics::TestUtils::randomInt64(1, 30 + 1 - 1);
            double dutyCycle = double(OpenMagnetics::TestUtils::randomInt64(1, 99 + 1 - 1)) / 100;
            double dcCurrent = 0;
            if (numberTurns.size() == 1) {
                dcCurrent = OpenMagnetics::TestUtils::randomInt(0, 30 - 1);
            }

            int waveShapeIndex = OpenMagnetics::TestUtils::randomInt(0, magic_enum::enum_count<WaveformLabel>() - 1);

            WaveformLabel waveShape = magic_enum::enum_cast<WaveformLabel>(waveShapeIndex).value();
            if (waveShape == WaveformLabel::BIPOLAR_RECTANGULAR ||
                waveShape == WaveformLabel::CUSTOM ||
                waveShape == WaveformLabel::RECTANGULAR ||
                waveShape == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME ||
                waveShape == WaveformLabel::UNIPOLAR_RECTANGULAR) {
                waveShape = WaveformLabel::TRIANGULAR;
            }

                for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
                    std::cout << "numberTurns: " << numberTurns[windingIndex] << std::endl;
                }
                for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
                    std::cout << "isolationSides: " << magic_enum::enum_name(isolationSides[windingIndex]) << std::endl;
                }
                std::cout << "frequency: " << frequency << std::endl;
                std::cout << "peakToPeak: " << peakToPeak << std::endl;
                std::cout << "magnetizingInductance: " << magnetizingInductance << std::endl;
                std::cout << "dutyCycle: " << dutyCycle << std::endl;
                std::cout << "dcCurrent: " << dcCurrent << std::endl;
                std::cout << "waveShape: " << magic_enum::enum_name(waveShape) << std::endl;

            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  waveShape,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  dcCurrent,
                                                                                                  turnsRatios);

            inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
            inputs.process();

            try {
                MagneticAdviser MagneticAdviser;
                auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 1);
                count--;
            }
            catch (...) {
                for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
                    std::cout << "numberTurns: " << numberTurns[windingIndex] << std::endl;
                }
                for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
                    std::cout << "isolationSides: " << magic_enum::enum_name(isolationSides[windingIndex]) << std::endl;
                }
                std::cout << "frequency: " << frequency << std::endl;
                std::cout << "peakToPeak: " << peakToPeak << std::endl;
                std::cout << "magnetizingInductance: " << magnetizingInductance << std::endl;
                std::cout << "dutyCycle: " << dutyCycle << std::endl;
                std::cout << "dcCurrent: " << dcCurrent << std::endl;
                std::cout << "waveShape: " << magic_enum::enum_name(waveShape) << std::endl;
                REQUIRE(false);
                return;
            }

        }
    }

    TEST_CASE("Test_MagneticAdviser_Random_0", "[adviser][magnetic-adviser][bug]") {

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {40, 75};
        std::vector<IsolationSide> isolationSides = {IsolationSide::SENARY,
                                                                    IsolationSide::SECONDARY};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 25280;
        double peakToPeak = 28;
        double magnetizingInductance = 0.004539;
        double dutyCycle = 0.68;
        double dcCurrent = 0;
        double temperature = 25;
        WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process();

        MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 1);

        if (masMagnetics.size() > 0) {
            auto masMagnetic = masMagnetics[0].first;
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_MagneticAdviser_Random_0_" + std::to_string(OpenMagnetics::TestUtils::randomInt(0, RAND_MAX)) + ".svg";
            outFile.append(filename);
            Painter painter(outFile, true);

            painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
            painter.paint_core(masMagnetic.get_mutable_magnetic());
            painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            #ifdef ENABLE_MATPLOTPP
            painter.export_svg();
            #else
                INFO("matplotplusplus disabled — skipping AdvancedPainter call");
            #endif

        }
    }

    TEST_CASE("Test_MagneticAdviser_Random_1", "[adviser][magnetic-adviser][bug]") {

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {24};
        std::vector<IsolationSide> isolationSides = {IsolationSide::SECONDARY};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 47405;
        double peakToPeak = 25;
        double magnetizingInductance = 0.000831;
        double dutyCycle = 0.05;
        double dcCurrent = 6;
        double temperature = 25;
        WaveformLabel waveShape = WaveformLabel::FLYBACK_SECONDARY;

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process();

        MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 1);

        if (masMagnetics.size() > 0) {
            auto masMagnetic = masMagnetics[0].first;
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_MagneticAdviser_Random_1_" + std::to_string(OpenMagnetics::TestUtils::randomInt(0, RAND_MAX)) + ".svg";
            outFile.append(filename);
            Painter painter(outFile, true);

            painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
            painter.paint_core(masMagnetic.get_mutable_magnetic());
            painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            #ifdef ENABLE_MATPLOTPP
            painter.export_svg();
            #else
                INFO("matplotplusplus disabled — skipping AdvancedPainter call");
            #endif

        }
    }

    TEST_CASE("Test_MagneticAdviser_Random_2", "[adviser][magnetic-adviser][bug]") {

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {45, 94};
        std::vector<IsolationSide> isolationSides = {IsolationSide::SECONDARY, IsolationSide::DENARY};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 569910;
        double peakToPeak = 1;
        double magnetizingInductance = 0.00354;
        double dutyCycle = 0.01;
        double dcCurrent = 0;
        double temperature = 25;
        WaveformLabel waveShape = WaveformLabel::TRIANGULAR;

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process();

        MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 1);

        if (masMagnetics.size() > 0) {
            auto masMagnetic = masMagnetics[0].first;
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_MagneticAdviser_Random_2_" + std::to_string(OpenMagnetics::TestUtils::randomInt(0, RAND_MAX)) + ".svg";
            outFile.append(filename);
            Painter painter(outFile, true);

            painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
            painter.paint_core(masMagnetic.get_mutable_magnetic());
            painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            #ifdef ENABLE_MATPLOTPP
            painter.export_svg();
            #else
                INFO("matplotplusplus disabled — skipping AdvancedPainter call");
            #endif

        }
    }

    TEST_CASE("Test_MagneticAdviser_Random_3", "[adviser][magnetic-adviser][bug]") {
        clear_databases();
        settings.reset();

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {78, 47, 1, 100};
        std::vector<IsolationSide> isolationSides = {IsolationSide::QUATERNARY,
                                                                    IsolationSide::TERTIARY,
                                                                    IsolationSide::QUINARY,
                                                                    IsolationSide::UNDENARY};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 811888;
        double peakToPeak = 1;
        double magnetizingInductance = 0.001213;
        double dutyCycle = 0.92;
        double dcCurrent = 0;
        double temperature = 25;
        WaveformLabel waveShape = WaveformLabel::FLYBACK_SECONDARY;

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process();

        MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 1);

        if (masMagnetics.size() > 0) {
            auto masMagnetic = masMagnetics[0].first;
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_MagneticAdviser_Random_3_" + std::to_string(OpenMagnetics::TestUtils::randomInt(0, RAND_MAX)) + ".svg";
            outFile.append(filename);
            Painter painter(outFile, true);

            painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
            painter.paint_core(masMagnetic.get_mutable_magnetic());
            painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            #ifdef ENABLE_MATPLOTPP
            painter.export_svg();
            #else
                INFO("matplotplusplus disabled — skipping AdvancedPainter call");
            #endif

        }
    }

    TEST_CASE("Test_MagneticAdviser_Web_0", "[adviser][magnetic-adviser][bug]") {
        settings.set_use_only_cores_in_stock(true);

        auto json_path_844 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_magneticadviser_web_0_844.json");
        std::ifstream json_file_844(json_path_844);
        OpenMagnetics::Inputs inputs = json::parse(json_file_844);

        MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 2);

        for (auto [mas, scoring] : masMagnetics) {
            std::string name = mas.get_magnetic().get_manufacturer_info().value().get_reference().value();
            auto masMagnetic = mas;
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            name = std::filesystem::path(std::regex_replace(std::string(name), std::regex(" "), "_")).string();
            name = std::filesystem::path(std::regex_replace(std::string(name), std::regex("\\,"), "_")).string();
            name = std::filesystem::path(std::regex_replace(std::string(name), std::regex("\\."), "_")).string();
            name = std::filesystem::path(std::regex_replace(std::string(name), std::regex("\\:"), "_")).string();
            name = std::filesystem::path(std::regex_replace(std::string(name), std::regex("\\/"), "_")).string();
            std::string filename = "Test_MagneticAdviser_Web_0_" + name + ".svg";
            outFile.append(filename);
            Painter painter(outFile, true);

            painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
            painter.paint_core(masMagnetic.get_mutable_magnetic());
            painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            #ifdef ENABLE_MATPLOTPP
            painter.export_svg();
            #else
                INFO("matplotplusplus disabled — skipping AdvancedPainter call");
            #endif

        }
    }

    TEST_CASE("Test_MagneticAdviser_Web_1", "[adviser][magnetic-adviser][bug]") {
        settings.set_use_only_cores_in_stock(true);

        OpenMagnetics::Inputs inputs = json::parse(R"({"designRequirements": {"insulation": {"altitude": {"maximum": 2000 }, "cti": "Group I", "pollutionDegree": "P2", "overvoltageCategory": "OVC-II", "insulationType": "Reinforced", "mainSupplyVoltage": {"maximum": 400 }, "standards": ["IEC 62368-1" ] }, "magnetizingInductance": {"nominal": 0.00039999999999999996 }, "name": "My Design Requirements", "turnsRatios": [{"nominal": 10 } ] }, "operatingPoints": [{"name": "Operating Point No. 1", "conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"ancillaryLabel": null, "data": [-2.5, 2.5, -2.5 ], "numberPeriods": null, "time": [0, 0.0000025, 0.00001 ] }, "processed": {"dutyCycle": 0.25, "peakToPeak": 5, "offset": 0, "label": "Triangular", "acEffectiveFrequency": 128062.87364573497, "effectiveFrequency": 128062.87364573497, "peak": 2.499999999999999, "rms": 1.4438454453780356, "thd": 0.37655984422983674 }, "harmonics": {"amplitudes": [5.289171878253285e-15, 1.9109142370405827, 0.6760173538930697, 0.21266521973828106, 1.0225829328048648e-15, 0.07680601066227859, 0.07559762456783126, 0.039376324985334776, 4.573645088368781e-16, 0.023974199977867476, 0.027568116438312952, 0.01617879664125947, 2.971250334652087e-16, 0.011696483630665971, 0.014340785210409782, 0.008885504288828753, 2.2404804069938507e-16, 0.0070081980332356445, 0.008903579153043517, 0.0056932375141851765, 2.0060103614117677e-16, 0.004737174081901557, 0.006158131362616469, 0.004020963170880106, 1.6456472022598636e-16, 0.0034711878538580374, 0.004586639500385065, 0.0030405807424970697, 8.639502758732636e-17, 0.0026976000170703556, 0.003608946997642845, 0.00242055067908846, 1.0574701291099621e-16, 0.0021941195434799954, 0.0029646243339663376, 0.00200725448389949, 1.078627841406342e-16, 0.001851823483504021, 0.002522860567684626, 0.0017217423565853373, 8.393610561532662e-17, 0.0016123930328628207, 0.002212325372118857, 0.0015202270093639689, 1.1562789535401302e-16, 0.0014424897474069642, 0.0019916921862528536, 0.0013770258476373952, 7.108157805209173e-17, 0.0013221392910760783, 0.0018359780830164152, 0.0012764919436558824, 9.44838447818231e-17, 0.0012390291425152763, 0.0017297261022555157, 0.0012089248059259444, 1.6792517082756795e-16, 0.0011855408939470595, 0.001663417261065822, 0.0011683976500078048, 8.640837577205282e-17, 0.0011571521557507119, 0.0016315323036584783, 0.0011515835097210575 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-533.3333333333333, 1599.9999999999998, 1599.9999999999998, -533.3333333333333, -533.3333333333333 ], "numberPeriods": null, "time": [0, 0, 0.0000025, 0.0000025, 0.00001 ] }, "processed": {"acEffectiveFrequency": 690362.0785055002, "average": -3.907985046680551e-13, "dutyCycle": 0.25, "effectiveFrequency": 690304.6703245766, "label": "Rectangular", "offset": 533.3333333333334, "peak": 1599.9999999999998, "peakToPeak": 2133.333333333333, "rms": 914.0872800534738, "thd": 0.9507084995254543 }, "harmonics": {"amplitudes": [16.666666666666664, 936.5743366559599, 678.5155874995734, 343.10398442523496, 33.33333333333333, 167.53211803338075, 224.71508018049997, 159.4087857982368, 33.33333333333333, 81.3927569897301, 133.0741261256697, 108.74207884995286, 33.33333333333333, 47.777866517707096, 93.16042574968293, 84.67373936814097, 33.33333333333333, 29.60335188747621, 70.47741191828814, 70.39622283854739, 33.33333333333333, 18.037041704179572, 55.61330685278376, 60.79168042183305, 33.33333333333333, 9.896923507307392, 44.944797116224166, 53.772679523931075, 33.33333333333333, 3.7542162789560147, 36.77766585778291, 48.326808642191196, 33.33333333333333, 1.1295053442673506, 30.211572300638238, 43.9020428459011, 33.33333333333333, 5.175841098510636, 24.721684875734596, 40.17025212715429, 33.33333333333333, 8.644545010278817, 19.979231122730887, 36.92259416373016, 33.33333333333333, 11.70597127653338, 15.765492529710798, 34.018188678401664, 33.33333333333333, 14.478186559464497, 11.926857377150895, 31.356776956938962, 33.33333333333333, 17.047466028773496, 8.349565339710303, 28.863097293035203, 33.33333333333333, 19.480403852519398, 4.944532917945068, 26.477336144614267, 33.33333333333333, 21.831581714888173, 1.6375616589827473, 24.148842611619322 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] } } }, {"current": {"harmonics": {"amplitudes": [0.0000010000000539855947, 19.037231223396237, 6.897984052452256, 2.326563601915805, 0.16964245489488294, 0.6827105539690119, 0.7639584293115499, 0.46512617224180786, 0.08523164109549347, 0.1853099567245175, 0.27319994312359214, 0.2033816061551562, 0.057281247031537384, 0.07597930993201436, 0.13794353928011985, 0.11736286786731481, 0.04345071383027455, 0.036440556746865825, 0.08222723635075783, 0.07810340034727467, 0.03527361775824142, 0.018254330742417037, 0.053961824922667094, 0.056682647846879576, 0.02992937223599544, 0.008573935375330551, 0.03763658083416089, 0.0435994370579927, 0.026210662569045903, 0.002867140660712981, 0.027318562880531572, 0.03496118462997451, 0.023515356873512087, 0.0007779637595615999, 0.020339562018896258, 0.02892190405184444, 0.021510534863982397, 0.0032750781058732345, 0.015353525869835107, 0.024509224663769844, 0.01999816717442588, 0.005104092422150859, 0.011619480338521122, 0.02116836728711663, 0.01885414149705068, 0.006537395429160153, 0.008699656143415562, 0.018561931478014347, 0.017997874963290205, 0.007740940114641076, 0.006318925086584056, 0.016473329652684467, 0.017376076319907764, 0.008823512955740787, 0.004293742261433124, 0.014756505062970141, 0.016953627505956508, 0.00986273033802779, 0.002493499634464789, 0.013307931595332303, 0.016708323451549476, 0.010919968097470808, 0.0008178599364280956, 0.01204994245895335 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] }, "processed": {"acEffectiveFrequency": 129404.9130926142, "average": -5.3734794391857577e-14, "dutyCycle": 0.2421874999999996, "effectiveFrequency": 129404.91309261405, "label": "Triangular", "offset": -0.000001, "peak": 24.999998999999956, "peakToPeak": 49.999999999999986, "rms": 14.438555983981303, "thd": 0.3878858672873039 }, "waveform": {"ancillaryLabel": null, "data": [-25.000000999999994, 24.999998999999992, -25.000000999999994 ], "numberPeriods": null, "time": [0, 0.000002421874999999996, 0.00001 ] } }, "frequency": 100000, "magneticFieldStrength": null, "magneticFluxDensity": null, "magnetizingCurrent": null, "name": "Primary winding excitation", "voltage": {"harmonics": {"amplitudes": [9.16666666666667, 114.01883606805356, 84.50783707683871, 45.4933446767725, 8.293206055601647, 17.4992480352357, 27.17398286677801, 22.031939878495663, 8.173210670026926, 6.406662607033046, 15.123340250417916, 15.15411908517724, 7.974502797768419, 1.9490359117489253, 9.560623024349578, 11.563504110131417, 7.698996104260727, 0.504510858012843, 6.182373666688366, 9.166959370632709, 7.349343869569635, 2.0540091893695203, 3.8205457805728327, 7.332930326725261, 6.928913435854551, 3.093876715391758, 2.0304195076924505, 5.808981955821909, 6.441753778022821, 3.8000361805623823, 0.6081443012940038, 4.4773238922780285, 5.892556509887897, 4.262766710863713, 0.5511898658329515, 3.2779337521433516, 5.286610701363717, 4.53337401311496, 1.5058617370415037, 2.178953825099018, 4.629751941830022, 4.643526917515465, 2.2899470149511534, 1.1635818333216585, 3.928306140216652, 4.614259582955264, 2.924044975741664, 0.22363428278401298, 3.189028603042417, 4.46053575052874, 3.420845617443641, 0.643874494083472, 2.4190389771205214, 4.1937095470147, 3.788199527266002, 1.4386641980705313, 1.6257526834677354, 3.822904217686624, 4.03087958389365, 2.1583263863852635, 0.8168095027463367, 3.3557907081369653, 4.151603816416468, 2.7990053199341562 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] }, "processed": {"acEffectiveFrequency": 707644.9123862737, "average": -1.6666666666666643, "dutyCycle": 0.2421874999999996, "effectiveFrequency": 706482.814815305, "label": "Unipolar Rectangular", "offset": -53.33333333333333, "peak": 213.33333333333331, "peakToPeak": 266.66666666666663, "rms": 113.33333333333348, "thd": 0.9813747916705575 }, "waveform": {"ancillaryLabel": null, "data": [-53.33333333333333, 213.33333333333331, 213.33333333333331, -53.33333333333333, -53.33333333333333 ], "numberPeriods": null, "time": [0, 0, 0.000002421874999999996, 0.000002421874999999996, 0.00001 ] } } } ] } ] })");
        MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 6);
        REQUIRE(masMagnetics.size() > 1);
    }

    TEST_CASE("Test_MagneticAdviser_Web_2", "[adviser][magnetic-adviser][bug]") {
        settings.set_use_only_cores_in_stock(false);

        OpenMagnetics::Inputs inputs = json::parse(R"({"designRequirements":{"magnetizingInductance":{"nominal":0.00001},"name":"My Design Requirements","turnsRatios":[]},"operatingPoints":[{"name":"Operating Point No. 1","conditions":{"ambientTemperature":42},"excitationsPerWinding":[{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.0000025,0.00001]},"processed":{"dutyCycle":0.25,"peakToPeak":10,"offset":0,"label":"Triangular","acEffectiveFrequency":128062.87364573497,"effectiveFrequency":128062.87364573497,"peak":4.999999999999998,"rms":2.887690890756071,"thd":0.37655984422983674},"harmonics":{"amplitudes":[1.057834375650657e-14,3.8218284740811654,1.3520347077861394,0.4253304394765621,2.0451658656097295e-15,0.15361202132455717,0.1511952491356625,0.07875264997066955,9.147290176737562e-16,0.04794839995573495,0.055136232876625904,0.03235759328251894,5.942500669304174e-16,0.023392967261331943,0.028681570420819563,0.017771008577657506,4.480960813987701e-16,0.014016396066471289,0.017807158306087034,0.011386475028370353,4.0120207228235354e-16,0.009474348163803114,0.012316262725232938,0.008041926341760212,3.291294404519727e-16,0.006942375707716075,0.00917327900077013,0.006081161484994139,1.7279005517465273e-16,0.005395200034140711,0.00721789399528569,0.00484110135817692,2.1149402582199242e-16,0.004388239086959991,0.005929248667932675,0.00401450896779898,2.157255682812684e-16,0.003703646967008042,0.005045721135369252,0.0034434847131706746,1.6787221123065324e-16,0.0032247860657256414,0.004424650744237714,0.0030404540187279378,2.3125579070802604e-16,0.0028849794948139283,0.003983384372505707,0.0027540516952747904,1.4216315610418345e-16,0.0026442785821521567,0.0036719561660328304,0.002552983887311765,1.889676895636462e-16,0.0024780582850305525,0.0034594522045110314,0.0024178496118518887,3.358503416551359e-16,0.002371081787894119,0.003326834522131644,0.0023367953000156096,1.7281675154410565e-16,0.0023143043115014237,0.0032630646073169567,0.002303167019442115],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}},"voltage":{"waveform":{"data":[-2.5,7.5,7.5,-2.5,-2.5],"time":[0,0,0.0000025,0.0000025,0.00001]},"processed":{"dutyCycle":0.25,"peakToPeak":10,"offset":0,"label":"Rectangular","acEffectiveFrequency":690362.0785055,"effectiveFrequency":690304.6703245763,"peak":7.5,"rms":4.284784125250653,"thd":0.9507084995254541},"harmonics":{"amplitudes":[0.078125,4.390192203074813,3.1805418164042503,1.6082999269932894,0.15625,0.7853068032814723,1.0533519383460939,0.747228683429235,0.15625,0.38152854838936,0.6237849662140766,0.5097284946091541,0.15625,0.22395874930175205,0.4366894957016387,0.3969081532881607,0.15625,0.1387657119725448,0.3303628683669758,0.32998229455569084,0.15625,0.08454863298834173,0.2606873758724238,0.2849610019773425,0.15625,0.046391828940503435,0.21067873648230084,0.25205943526842706,0.15625,0.017597888807606342,0.1723953087083574,0.22653191551027138,0.15625,0.005294556301253328,0.14161674515924202,0.20579082584016142,0.15625,0.02426175514926864,0.11588289785500591,0.18829805684603593,0.15625,0.040521304735682134,0.09365264588780098,0.17307466014248524,0.15625,0.05487174035875039,0.07390074623301954,0.1594602594300079,0.15625,0.06786649949748992,0.05590714395539467,0.1469848919856514,0.15625,0.07990999700987583,0.03913858752989208,0.13529576856110245,0.15625,0.09131439305868466,0.023177498052867507,0.1241125131778793,0.15625,0.10233553928853857,0.007676070276481761,0.11319769974196572],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}}}]}]})");
        MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 6);
        REQUIRE(masMagnetics.size() > 1);
    }

    TEST_CASE("Test_MagneticAdviser_Web_3", "[adviser][magnetic-adviser][bug]") {
        settings.set_use_only_cores_in_stock(true);

        OpenMagnetics::Inputs inputs = json::parse(R"({"designRequirements": {"isolationSides": ["primary" ], "magnetizingInductance": {"nominal": 0.000049999999999999996 }, "turnsRatios": [], "name": "My Design Requirements" }, "operatingPoints": [{"conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"current": null, "frequency": 100000, "magneticFluxDensity": {"processed": {"average": -0.0008402065119815163, "dutyCycle": 0.25, "label": "Custom", "offset": 0.00025834764404089605, "peak": 0.11157850212098047, "peakToPeak": 0.22264030895387915 }, "waveform": {"data": [-0.10403106023435499, -0.10637464243386952, 0.11157850212098047, -0.11106180683289868, -0.10403106023435499 ], "time": [0, 7.8125e-8, 0.0000025000000000000015, 0.000009921875000000021, 0.000010000000000000021 ] } }, "magnetizingCurrent": {"harmonics": {"amplitudes": [0.0004152223659327984, 0.042441838144660536, 0.01496020903582336, 0.0046781309589203806, 0.00007238605418457529, 0.0016575407471086866, 0.00160867800926016, 0.0008260138819145056, 0.00007238605418458149, 0.00048463999871970116, 0.0005406460365225253, 0.0003127993597737691, 0.00007238605418458021, 0.00021508894336653934, 0.0002465099269605511, 0.00015509187555504158, 0.00007238605418458009, 0.00011650268024990807, 0.00012560279845175753, 0.00009114503667276431, 0.00007238605418458025, 0.00007451661945210314, 0.00006455228734198571, 0.00006388554543154429, 0.00007238605418458002, 0.00005741153788913295, 0.0000296070293303051, 0.000053756608745154496, 0.00007238605418458025, 0.000051935950882856886, 0.000007866085208646598, 0.00005125276440434498, 0.00007238605418458002, 0.00005124062704727633, 0.000006461715186690873, 0.000051601981925933585, 0.00007238605418457974, 0.000052153466095872496, 0.000016285214222190135, 0.00005278481183552113, 0.00007238605418458, 0.000053431039924951465, 0.000023190584075869416, 0.000054054727364434106, 0.00007238605418457993, 0.00005463503709822476, 0.00002809680331335256, 0.00005516102161109927, 0.00007238605418457995, 0.00005562755227300169, 0.000031559417180106934, 0.00005603284378578814, 0.0000723860541845798, 0.00005637694553791417, 0.000033922142034865184, 0.0000566608204897135, 0.00007238605418457978, 0.000056885782643672325, 0.000035396651460918886, 0.000057053154515285994, 0.00007238605418458014, 0.00005716406046156705, 0.00003610567712631751, 0.00005721930487623161 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] }, "processed": {"average": -0.0004152223659328003, "dutyCycle": 0.25, "label": "Custom", "offset": 0.00012767304045150502, "peak": 0.0551410742207318, "peakToPeak": 0.11002680236056059 }, "waveform": {"data": [-0.05141119753896891, -0.05256937440592218, -0.04909484380506238, -0.04562031320420257, -0.042145782603342764, -0.03867125200248296, -0.03519672140162316, -0.031722190800763356, -0.028247660199903553, -0.024773129599043753, -0.02129859899818395, -0.017824068397324148, -0.014349537796464347, -0.010875007195604539, -0.007400476594744737, -0.00392594599388493, -0.0004514153930251279, 0.003023115207834674, 0.006497645808694477, 0.009972176409554279, 0.01344670701041408, 0.01692123761127388, 0.020395768212133684, 0.023870298812993487, 0.02734482941385329, 0.030819360014713092, 0.034293890615572906, 0.037768421216432715, 0.041242951817292525, 0.04471748241815233, 0.04819201301901215, 0.051666543619871975, 0.0551410742207318, 0.053982897353778536, 0.05282472048682527, 0.05166654361987202, 0.050508366752918754, 0.04935018988596549, 0.048192013019012235, 0.04703383615205897, 0.04587565928510571, 0.04471748241815245, 0.04355930555119919, 0.04240112868424593, 0.04124295181729267, 0.04008477495033941, 0.03892659808338614, 0.03776842121643287, 0.0366102443494796, 0.03545206748252632, 0.03429389061557305, 0.03313571374861978, 0.03197753688166651, 0.03081936001471324, 0.02966118314775997, 0.0285030062808067, 0.02734482941385343, 0.02618665254690016, 0.02502847567994689, 0.02387029881299362, 0.022712121946040345, 0.021553945079087075, 0.020395768212133805, 0.019237591345180535, 0.018079414478227265, 0.016921237611273992, 0.015763060744320722, 0.014604883877367452, 0.01344670701041418, 0.01228853014346091, 0.011130353276507639, 0.009972176409554369, 0.008813999542601097, 0.007655822675647827, 0.0064976458086945564, 0.005339468941741286, 0.004181292074788015, 0.0030231152078347445, 0.0018649383408814739, 0.0007067614739282031, -0.00045141539302506766, -0.0016095922599783383, -0.002767769126931609, -0.00392594599388488, -0.005084122860838151, -0.0062422997277914215, -0.007400476594744692, -0.008558653461697962, -0.009716830328651234, -0.010875007195604504, -0.012033184062557776, -0.013191360929511045, -0.014349537796464315, -0.015507714663417587, -0.016665891530370857, -0.017824068397324127, -0.0189822452642774, -0.02014042213123067, -0.02129859899818394, -0.02245677586513721, -0.02361495273209048, -0.024773129599043753, -0.025931306465997023, -0.027089483332950293, -0.028247660199903563, -0.029405837066856837, -0.030564013933810107, -0.03172219080076338, -0.032880367667716646, -0.034038544534669916, -0.035196721401623186, -0.036354898268576456, -0.03751307513552973, -0.038671252002483, -0.03982942886943627, -0.04098760573638954, -0.04214578260334281, -0.04330395947029608, -0.04446213633724935, -0.04562031320420262, -0.04677849007115589, -0.04793666693810917, -0.04909484380506244, -0.05025302067201571, -0.05141119753896898, -0.05256937440592225, -0.05372755127287552, -0.05488572813982879 ], "time": [0, 7.8125e-8, 1.5625e-7, 2.3437500000000003e-7, 3.125e-7, 3.90625e-7, 4.6875e-7, 5.468750000000001e-7, 6.25e-7, 7.03125e-7, 7.8125e-7, 8.59375e-7, 9.375e-7, 0.0000010156250000000001, 0.0000010937500000000001, 0.0000011718750000000001, 0.00000125, 0.000001328125, 0.00000140625, 0.000001484375, 0.0000015625, 0.000001640625, 0.00000171875, 0.000001796875, 0.000001875, 0.000001953125, 0.0000020312500000000002, 0.0000021093750000000005, 0.0000021875000000000007, 0.000002265625000000001, 0.000002343750000000001, 0.0000024218750000000013, 0.0000025000000000000015, 0.0000025781250000000017, 0.000002656250000000002, 0.000002734375000000002, 0.0000028125000000000023, 0.0000028906250000000025, 0.0000029687500000000027, 0.000003046875000000003, 0.000003125000000000003, 0.0000032031250000000033, 0.0000032812500000000035, 0.0000033593750000000037, 0.000003437500000000004, 0.000003515625000000004, 0.0000035937500000000043, 0.0000036718750000000045, 0.0000037500000000000048, 0.000003828125000000005, 0.000003906250000000005, 0.000003984375000000005, 0.0000040625000000000056, 0.000004140625000000006, 0.000004218750000000006, 0.000004296875000000006, 0.000004375000000000006, 0.000004453125000000007, 0.000004531250000000007, 0.000004609375000000007, 0.000004687500000000007, 0.000004765625000000007, 0.000004843750000000008, 0.000004921875000000008, 0.000005000000000000008, 0.000005078125000000008, 0.0000051562500000000084, 0.000005234375000000009, 0.000005312500000000009, 0.000005390625000000009, 0.000005468750000000009, 0.0000055468750000000095, 0.00000562500000000001, 0.00000570312500000001, 0.00000578125000000001, 0.00000585937500000001, 0.0000059375000000000105, 0.000006015625000000011, 0.000006093750000000011, 0.000006171875000000011, 0.000006250000000000011, 0.0000063281250000000115, 0.000006406250000000012, 0.000006484375000000012, 0.000006562500000000012, 0.000006640625000000012, 0.0000067187500000000125, 0.000006796875000000013, 0.000006875000000000013, 0.000006953125000000013, 0.000007031250000000013, 0.0000071093750000000136, 0.000007187500000000014, 0.000007265625000000014, 0.000007343750000000014, 0.000007421875000000014, 0.000007500000000000015, 0.000007578125000000015, 0.000007656250000000015, 0.000007734375000000015, 0.000007812500000000015, 0.000007890625000000016, 0.000007968750000000016, 0.000008046875000000016, 0.000008125000000000016, 0.000008203125000000016, 0.000008281250000000017, 0.000008359375000000017, 0.000008437500000000017, 0.000008515625000000017, 0.000008593750000000017, 0.000008671875000000018, 0.000008750000000000018, 0.000008828125000000018, 0.000008906250000000018, 0.000008984375000000018, 0.000009062500000000019, 0.000009140625000000019, 0.000009218750000000019, 0.00000929687500000002, 0.00000937500000000002, 0.00000945312500000002, 0.00000953125000000002, 0.00000960937500000002, 0.00000968750000000002, 0.00000976562500000002, 0.00000984375000000002, 0.000009921875000000021 ] } }, "name": "Primary winding excitation", "voltage": {"harmonics": {"amplitudes": [0.078125, 4.390192203074813, 3.1805418164042503, 1.6082999269932894, 0.15625, 0.7853068032814723, 1.0533519383460939, 0.747228683429235, 0.15625, 0.38152854838936, 0.6237849662140766, 0.5097284946091541, 0.15625, 0.22395874930175205, 0.4366894957016387, 0.3969081532881607, 0.15625, 0.1387657119725448, 0.3303628683669758, 0.32998229455569084, 0.15625, 0.08454863298834173, 0.2606873758724238, 0.2849610019773425, 0.15625, 0.046391828940503435, 0.21067873648230084, 0.25205943526842706, 0.15625, 0.017597888807606342, 0.1723953087083574, 0.22653191551027138, 0.15625, 0.005294556301253328, 0.14161674515924202, 0.20579082584016142, 0.15625, 0.02426175514926864, 0.11588289785500591, 0.18829805684603593, 0.15625, 0.040521304735682134, 0.09365264588780098, 0.17307466014248524, 0.15625, 0.05487174035875039, 0.07390074623301954, 0.1594602594300079, 0.15625, 0.06786649949748992, 0.05590714395539467, 0.1469848919856514, 0.15625, 0.07990999700987583, 0.03913858752989208, 0.13529576856110245, 0.15625, 0.09131439305868466, 0.023177498052867507, 0.1241125131778793, 0.15625, 0.10233553928853857, 0.007676070276481761, 0.11319769974196572 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] }, "processed": {"acEffectiveFrequency": 690362.0785055, "average": -0.078125, "dutyCycle": 0.2421874999999996, "effectiveFrequency": 690304.6703245763, "label": "Rectangular", "offset": 0, "peak": 7.5, "peakToPeak": 10, "rms": 4.284784125250653, "thd": 0.9507084995254541 }, "waveform": {"data": [-2.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5 ], "time": [0, 7.8125e-8, 1.5625e-7, 2.3437500000000003e-7, 3.125e-7, 3.90625e-7, 4.6875e-7, 5.468750000000001e-7, 6.25e-7, 7.03125e-7, 7.8125e-7, 8.59375e-7, 9.375e-7, 0.0000010156250000000001, 0.0000010937500000000001, 0.0000011718750000000001, 0.00000125, 0.000001328125, 0.00000140625, 0.000001484375, 0.0000015625, 0.000001640625, 0.00000171875, 0.000001796875, 0.000001875, 0.000001953125, 0.0000020312500000000002, 0.0000021093750000000005, 0.0000021875000000000007, 0.000002265625000000001, 0.000002343750000000001, 0.0000024218750000000013, 0.0000025000000000000015, 0.0000025781250000000017, 0.000002656250000000002, 0.000002734375000000002, 0.0000028125000000000023, 0.0000028906250000000025, 0.0000029687500000000027, 0.000003046875000000003, 0.000003125000000000003, 0.0000032031250000000033, 0.0000032812500000000035, 0.0000033593750000000037, 0.000003437500000000004, 0.000003515625000000004, 0.0000035937500000000043, 0.0000036718750000000045, 0.0000037500000000000048, 0.000003828125000000005, 0.000003906250000000005, 0.000003984375000000005, 0.0000040625000000000056, 0.000004140625000000006, 0.000004218750000000006, 0.000004296875000000006, 0.000004375000000000006, 0.000004453125000000007, 0.000004531250000000007, 0.000004609375000000007, 0.000004687500000000007, 0.000004765625000000007, 0.000004843750000000008, 0.000004921875000000008, 0.000005000000000000008, 0.000005078125000000008, 0.0000051562500000000084, 0.000005234375000000009, 0.000005312500000000009, 0.000005390625000000009, 0.000005468750000000009, 0.0000055468750000000095, 0.00000562500000000001, 0.00000570312500000001, 0.00000578125000000001, 0.00000585937500000001, 0.0000059375000000000105, 0.000006015625000000011, 0.000006093750000000011, 0.000006171875000000011, 0.000006250000000000011, 0.0000063281250000000115, 0.000006406250000000012, 0.000006484375000000012, 0.000006562500000000012, 0.000006640625000000012, 0.0000067187500000000125, 0.000006796875000000013, 0.000006875000000000013, 0.000006953125000000013, 0.000007031250000000013, 0.0000071093750000000136, 0.000007187500000000014, 0.000007265625000000014, 0.000007343750000000014, 0.000007421875000000014, 0.000007500000000000015, 0.000007578125000000015, 0.000007656250000000015, 0.000007734375000000015, 0.000007812500000000015, 0.000007890625000000016, 0.000007968750000000016, 0.000008046875000000016, 0.000008125000000000016, 0.000008203125000000016, 0.000008281250000000017, 0.000008359375000000017, 0.000008437500000000017, 0.000008515625000000017, 0.000008593750000000017, 0.000008671875000000018, 0.000008750000000000018, 0.000008828125000000018, 0.000008906250000000018, 0.000008984375000000018, 0.000009062500000000019, 0.000009140625000000019, 0.000009218750000000019, 0.00000929687500000002, 0.00000937500000000002, 0.00000945312500000002, 0.00000953125000000002, 0.00000960937500000002, 0.00000968750000000002, 0.00000976562500000002, 0.00000984375000000002, 0.000009921875000000021 ] } } } ], "name": "Operating Point No. 1" } ]})");
        MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 6);

        REQUIRE(masMagnetics.size() > 1);
        for (auto [mas, scoring] : masMagnetics) {
            REQUIRE_THAT(mas.get_outputs()[0].get_inductance()->get_magnetizing_inductance().get_magnetizing_inductance().get_nominal().value(), Catch::Matchers::WithinAbs(mas.get_inputs().get_design_requirements().get_magnetizing_inductance().get_nominal().value(), mas.get_inputs().get_design_requirements().get_magnetizing_inductance().get_nominal().value() * 0.25));
        }
    }

    TEST_CASE("Test_MagneticAdviser_Web_4", "[adviser][magnetic-adviser][bug]") {
        settings.set_use_only_cores_in_stock(true);

        auto json_path_906 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_magneticadviser_web_4_906.json");
        std::ifstream json_file_906(json_path_906);
        OpenMagnetics::Inputs inputs = json::parse(json_file_906);
        MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 6);

        REQUIRE(masMagnetics.size() > 1);
        for (auto [mas, scoring] : masMagnetics) {
            REQUIRE_THAT(mas.get_outputs()[0].get_inductance()->get_magnetizing_inductance().get_magnetizing_inductance().get_nominal().value(), Catch::Matchers::WithinAbs(mas.get_inputs().get_design_requirements().get_magnetizing_inductance().get_nominal().value(), mas.get_inputs().get_design_requirements().get_magnetizing_inductance().get_nominal().value() * 0.25));
            std::string name = mas.get_magnetic().get_core().get_name().value();
            auto masMagnetic = mas;
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
        }
    }

    TEST_CASE("Test_MagneticAdviser_Web_5", "[adviser][magnetic-adviser][bug]") {
        settings.set_use_only_cores_in_stock(true);

        auto json_path_922 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_magneticadviser_web_5_922.json");
        std::ifstream json_file_922(json_path_922);
        OpenMagnetics::Inputs inputs = json::parse(json_file_922);
        MagneticAdviser MagneticAdviser;
        std::map<MagneticFilters, double> weights;
        weights[MagneticFilters::COST] = 28.57142857142857;
        weights[MagneticFilters::LOSSES] = 21.428571428571427;
        weights[MagneticFilters::DIMENSIONS] = 50;

        
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, weights, 6);

        REQUIRE(masMagnetics.size() > 1);
        for (auto [mas, scoring] : masMagnetics) {
            REQUIRE_THAT(mas.get_outputs()[0].get_inductance()->get_magnetizing_inductance().get_magnetizing_inductance().get_nominal().value(), Catch::Matchers::WithinAbs(mas.get_inputs().get_design_requirements().get_magnetizing_inductance().get_nominal().value(), mas.get_inputs().get_design_requirements().get_magnetizing_inductance().get_nominal().value() * 0.25));
            std::string name = mas.get_magnetic().get_core().get_name().value();
            auto masMagnetic = mas;
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
        }
    }

    TEST_CASE("Test_MagneticAdviser_Web_6", "[adviser][magnetic-adviser][bug]") {
        settings.set_use_only_cores_in_stock(true);

        auto json_path_946 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_magneticadviser_web_6_946.json");
        std::ifstream json_file_946(json_path_946);
        OpenMagnetics::Inputs inputs = json::parse(json_file_946);
        MagneticAdviser MagneticAdviser;
        std::map<MagneticFilters, double> weights;
        weights[MagneticFilters::COST] = 28.57142857142857;
        weights[MagneticFilters::LOSSES] = 21.428571428571427;
        weights[MagneticFilters::DIMENSIONS] = 50;

        
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, weights, 6);

        REQUIRE(masMagnetics.size() > 1);
        for (auto [mas, scoring] : masMagnetics) {
            REQUIRE_THAT(mas.get_outputs()[0].get_inductance()->get_magnetizing_inductance().get_magnetizing_inductance().get_nominal().value(), Catch::Matchers::WithinAbs(mas.get_inputs().get_design_requirements().get_magnetizing_inductance().get_nominal().value(), mas.get_inputs().get_design_requirements().get_magnetizing_inductance().get_nominal().value() * 0.25));
            std::string name = mas.get_magnetic().get_core().get_name().value();
            auto masMagnetic = mas;
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
        }
    }

    TEST_CASE("Test_MagneticAdviser_Web_7", "[adviser][magnetic-adviser][bug]") {
        settings.set_use_only_cores_in_stock(true);

        OpenMagnetics::Inputs inputs = json::parse(R"({"designRequirements":{"magnetizingInductance":{"nominal": 0.0001},"minimumImpedance":[{"frequency":100000,"impedance":{"magnitude":1000}}],"name":"My Design Requirements","turnsRatios":[{"nominal":1}]},"operatingPoints":[{"name":"Operating Point No. 1","conditions":{"ambientTemperature":42},"excitationsPerWinding":[{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular","acEffectiveFrequency":110746.40291779551,"effectiveFrequency":110746.40291779551,"peak":5,"rms":2.8874560332150576,"thd":0.12151487440704967},"harmonics":{"amplitudes":[1.1608769501236793e-14,4.05366124583194,1.787369544444173e-15,0.4511310569983995,9.749053004706756e-16,0.16293015292554872,4.036157626725542e-16,0.08352979924600704,3.4998295008010614e-16,0.0508569581336163,3.1489164048780735e-16,0.034320410449418075,3.142469873118059e-16,0.024811988673843106,2.3653352035940994e-16,0.018849001010678823,2.9306524147249266e-16,0.014866633059596499,1.796485796132569e-16,0.012077180559557785,1.6247782523152451e-16,0.010049063750920609,1.5324769149805092e-16,0.008529750975091871,1.0558579733068502e-16,0.007363501410705499,7.513269775674661e-17,0.006450045785294609,5.871414177162291e-17,0.005722473794997712,9.294731722001391e-17,0.005134763398167541,1.194820309200107e-16,0.004654430423785411,8.2422739080512e-17,0.004258029771397032,9.5067306351894e-17,0.0039283108282380024,1.7540347128474968e-16,0.0036523670873925395,9.623794010508822e-17,0.0034204021424253787,1.4083470894369491e-16,0.0032248884817922927,1.4749333016985644e-16,0.0030599828465501895,1.0448590642474364e-16,0.002921112944200383,7.575487373767413e-18,0.002804680975178716,7.419510610361002e-17,0.0027078483284668376,3.924741709073613e-17,0.0026283777262804953,2.2684279102637236e-17,0.0025645167846443107,8.997077625295079e-17,0.0025149120164513483,7.131074184849219e-17,0.0024785457043284276,9.354417496250849e-17,0.0024546904085875065,1.2488589642405877e-16,0.0024428775264784264],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}},"voltage":{"waveform":{"data":[-20.5,70.5,70.5,-20.5,-20.5],"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular","acEffectiveFrequency":591485.5360118389,"effectiveFrequency":553357.3374711702,"peak":70.5,"rms":51.572309672924284,"thd":0.4833151484524849},"harmonics":{"amplitudes":[24.2890625,57.92076613061847,1.421875,19.27588907896988,1.421875,11.528257939603122,1.421875,8.194467538528329,1.421875,6.331896912839248,1.421875,5.137996046859012,1.421875,4.304077056139349,1.421875,3.6860723299088454,1.421875,3.207698601961777,1.421875,2.8247804703632298,1.421875,2.509960393415185,1.421875,2.2453859950684323,1.421875,2.01890737840567,1.421875,1.8219644341144872,1.421875,1.6483482744897402,1.421875,1.4934420157473332,1.421875,1.3537375367153817,1.421875,1.2265178099275544,1.421875,1.1096421410704556,1.421875,1.0013973584174929,1.421875,0.9003924136274832,1.421875,0.8054822382572133,1.421875,0.7157117294021269,1.421875,0.6302738400635857,1.421875,0.5484777114167545,1.421875,0.46972405216147894,1.421875,0.3934858059809043,1.421875,0.31929270856030145,1.421875,0.24671871675852053,1.421875,0.17537155450693565,1.421875,0.10488380107099537,1.421875,0.034905072061178544],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}}},{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"ancillaryLabel":null,"data":[-5,5,-5],"numberPeriods":null,"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular","acEffectiveFrequency":110746.40291779551,"effectiveFrequency":110746.40291779551,"peak":5,"rms":2.8874560332150576,"thd":0.12151487440704967},"harmonics":{"amplitudes":[1.1608769501236793e-14,4.05366124583194,1.787369544444173e-15,0.4511310569983995,9.749053004706756e-16,0.16293015292554872,4.036157626725542e-16,0.08352979924600704,3.4998295008010614e-16,0.0508569581336163,3.1489164048780735e-16,0.034320410449418075,3.142469873118059e-16,0.024811988673843106,2.3653352035940994e-16,0.018849001010678823,2.9306524147249266e-16,0.014866633059596499,1.796485796132569e-16,0.012077180559557785,1.6247782523152451e-16,0.010049063750920609,1.5324769149805092e-16,0.008529750975091871,1.0558579733068502e-16,0.007363501410705499,7.513269775674661e-17,0.006450045785294609,5.871414177162291e-17,0.005722473794997712,9.294731722001391e-17,0.005134763398167541,1.194820309200107e-16,0.004654430423785411,8.2422739080512e-17,0.004258029771397032,9.5067306351894e-17,0.0039283108282380024,1.7540347128474968e-16,0.0036523670873925395,9.623794010508822e-17,0.0034204021424253787,1.4083470894369491e-16,0.0032248884817922927,1.4749333016985644e-16,0.0030599828465501895,1.0448590642474364e-16,0.002921112944200383,7.575487373767413e-18,0.002804680975178716,7.419510610361002e-17,0.0027078483284668376,3.924741709073613e-17,0.0026283777262804953,2.2684279102637236e-17,0.0025645167846443107,8.997077625295079e-17,0.0025149120164513483,7.131074184849219e-17,0.0024785457043284276,9.354417496250849e-17,0.0024546904085875065,1.2488589642405877e-16,0.0024428775264784264],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}},"voltage":{"waveform":{"ancillaryLabel":null,"data":[-50,50,50,-50,-50],"numberPeriods":null,"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular","acEffectiveFrequency":591485.536011839,"effectiveFrequency":591449.4202715514,"peak":50,"rms":50,"thd":0.48331514845248497},"harmonics":{"amplitudes":[0.78125,63.64919355013018,1.5625,21.18229569117569,1.5625,12.668415318245188,1.5625,9.004909382998164,1.5625,6.958128475647527,1.5625,5.646149502042871,1.5625,4.729755006746538,1.5625,4.050628933965765,1.5625,3.524943518639316,1.5625,3.104154363036517,1.5625,2.7581982345221827,1.5625,2.467457137437843,1.5625,2.2185795367095267,1.5625,2.0021587188071255,1.5625,1.8113717302085082,1.5625,1.6411450722498175,1.5625,1.487623666720196,1.5625,1.3478217691511587,1.5625,1.2193869682092893,1.5625,1.100436657601639,1.5625,0.9894422127774558,1.5625,0.8851453167661671,1.5625,0.7864964059364037,1.5625,0.6926086154544899,1.5625,0.60272275979863,1.5625,0.5161802771005264,1.5625,0.43240198459440116,1.5625,0.3508711083080249,1.5625,0.27111946896540395,1.5625,0.192715993963664,1.5625,0.11525692425384548,1.5625,0.03835722204524927],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}}},{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular"}},"voltage":{"waveform":{"data":[-20.5,70.5,70.5,-20.5,-20.5],"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular"}}},{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular"}},"voltage":{"waveform":{"data":[-20.5,70.5,70.5,-20.5,-20.5],"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular"}}}]}]})");
        MagneticAdviser MagneticAdviser;
        std::map<MagneticFilters, double> weights;
        weights[MagneticFilters::COST] = 28.57142857142857;
        weights[MagneticFilters::LOSSES] = 21.428571428571427;
        weights[MagneticFilters::DIMENSIONS] = 50;

        
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, weights, 6);

        if (plot && masMagnetics.size() > 0) {
            auto masMagnetic = masMagnetics[0].first;
            // MagneticAdviser::preview_magnetic(masMagnetic);
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_MagneticAdviser_Web_7.svg";
            outFile.append(filename);
            Painter painter(outFile, true);

            painter.paint_core(masMagnetic.get_mutable_magnetic());
            painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            #ifdef ENABLE_MATPLOTPP
            painter.export_svg();
            #else
                INFO("matplotplusplus disabled — skipping AdvancedPainter call");
            #endif

        }
    }

    TEST_CASE("Test_MagneticAdviser_Inductor", "[adviser][magnetic-adviser]") {

        std::vector<double> turnsRatios;

        std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY};

        double frequency = 50000;
        double peakToPeak = 90;
        double magnetizingInductance = 0.000110;
        double dutyCycle = 0.5;
        double dcCurrent = 0;
        double temperature = 25;
        WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process();

        MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 1);

        if (plot && masMagnetics.size() > 0) {
            auto masMagnetic = masMagnetics[0].first;
            // MagneticAdviser::preview_magnetic(masMagnetic);
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_MagneticAdviser_Inductor.svg";
            outFile.append(filename);
            Painter painter(outFile, true);

            painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
            painter.paint_core(masMagnetic.get_mutable_magnetic());
            painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            #ifdef ENABLE_MATPLOTPP
            painter.export_svg();
            #else
                INFO("matplotplusplus disabled — skipping AdvancedPainter call");
            #endif

        }
    }

    TEST_CASE("Test_MagneticAdviser_Inductor_Only_Toroidal_Cores", "[adviser][magnetic-adviser][smoke-test]") {
        settings.reset();
        clear_databases();
        settings.set_use_toroidal_cores(true);
        settings.set_use_concentric_cores(false);


        std::vector<double> turnsRatios;

        std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY};

        double frequency = 50000;
        double peakToPeak = 90;
        double magnetizingInductance = 0.000110;
        double dutyCycle = 0.5;
        double dcCurrent = 0;
        double temperature = 25;
        WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process();

        MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 1);

        REQUIRE(masMagnetics.size() > 0);

        if (plot)
        for (auto [masMagnetic, scoring] : masMagnetics) {
            // MagneticAdviser::preview_magnetic(masMagnetic);
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_MagneticAdviser_Inductor_Only_Toroidal_Cores.svg";
            outFile.append(filename);
            Painter painter(outFile, true);

            painter.paint_core(masMagnetic.get_mutable_magnetic());
            // painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            #ifdef ENABLE_MATPLOTPP
            painter.export_svg();
            #else
                INFO("matplotplusplus disabled — skipping AdvancedPainter call");
            #endif

        }
    }

    TEST_CASE("Test_IEEE_Article_0", "[adviser][magnetic-adviser]") {
        settings.reset();
        clear_databases();
        settings.set_use_only_cores_in_stock(false);
        settings.set_wire_adviser_include_foil(true);
        settings.set_use_toroidal_cores(false);
        settings.set_use_concentric_cores(true);
        json inputsJson;

        inputsJson["operatingPoints"] = json::array();
        json operatingPointJson = json();
        operatingPointJson["name"] = "Nominal";
        operatingPointJson["conditions"] = json();
        operatingPointJson["conditions"]["ambientTemperature"] = 100;

        operatingPointJson["excitationsPerWinding"] = json::array();
        {
            json windingExcitation = json();
            windingExcitation["frequency"] = 115000;
            windingExcitation["current"]["processed"]["label"] = "Sinusoidal";
            windingExcitation["current"]["processed"]["peakToPeak"] = 4 * 2;
            windingExcitation["current"]["processed"]["offset"] = 0;
            windingExcitation["current"]["processed"]["dutyCycle"] = 0.5;
            windingExcitation["voltage"]["processed"]["label"] = "Rectangular";
            windingExcitation["voltage"]["processed"]["peakToPeak"] = 380;
            windingExcitation["voltage"]["processed"]["offset"] = 0;
            windingExcitation["voltage"]["processed"]["dutyCycle"] = 0.5;
            operatingPointJson["excitationsPerWinding"].push_back(windingExcitation);
        }
        {
            json windingExcitation = json();
            windingExcitation["frequency"] = 115000;
            windingExcitation["current"]["processed"]["label"] = "Sinusoidal";
            windingExcitation["current"]["processed"]["peakToPeak"] = 64 * 2;
            windingExcitation["current"]["processed"]["offset"] = 0;
            windingExcitation["current"]["processed"]["dutyCycle"] = 0.5;
            windingExcitation["voltage"]["processed"]["label"] = "Rectangular";
            windingExcitation["voltage"]["processed"]["peakToPeak"] = 24;
            windingExcitation["voltage"]["processed"]["offset"] = 0;
            windingExcitation["voltage"]["processed"]["dutyCycle"] = 0.5;
            operatingPointJson["excitationsPerWinding"].push_back(windingExcitation);
        }
        inputsJson["operatingPoints"].push_back(operatingPointJson);

        inputsJson["designRequirements"] = json();
        inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 250e-6;
        inputsJson["designRequirements"]["turnsRatios"] = json::array();
        json turnRatioSecondary;
        turnRatioSecondary["minimum"] = 16;
        turnRatioSecondary["maximum"] = 16;
        turnRatioSecondary["nominal"] = 16;
        inputsJson["designRequirements"]["turnsRatios"].push_back(turnRatioSecondary);


        OpenMagnetics::Inputs inputs(inputsJson);
        std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY,
                                                                    IsolationSide::SECONDARY};


        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process();
        MagneticAdviser MagneticAdviser;

        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 3);

        REQUIRE(masMagnetics.size() > 0);

        for (auto [masMagnetic, scoring] : masMagnetics) {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");

            {
                auto outFile = outputFilePath;
                std::string filename = "Test_IEEE_Article_0_" +  std::to_string(scoring) + ".mas.json";
                outFile.append(filename);
                to_file(outFile, masMagnetic);
            }

            // MagneticAdviser::preview_magnetic(masMagnetic);
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outFile = outputFilePath;
            std::string filename = "Test_IEEE_Article_0_" + std::to_string(scoring) + ".svg";
            outFile.append(filename);
            Painter painter(outFile, true);

            settings.set_painter_mode(PainterModes::CONTOUR);
            settings.set_painter_logarithmic_scale(true);
            settings.set_painter_include_fringing(true);
            settings.set_painter_maximum_value_colorbar(std::nullopt);
            settings.set_painter_minimum_value_colorbar(std::nullopt);
            // settings.set_painter_number_points_x(200);
            // settings.set_painter_number_points_y(200);
            painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
            painter.paint_core(masMagnetic.get_mutable_magnetic());
            // painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            #ifdef ENABLE_MATPLOTPP
            painter.export_svg();
            #else
                INFO("matplotplusplus disabled — skipping AdvancedPainter call");
            #endif

        }
    }

    TEST_CASE("Test_Magnetic_Adviser_Psim", "[adviser][magnetic-adviser][smoke-test]") {
        settings.reset();
        clear_databases();
        settings.set_use_only_cores_in_stock(true);
        settings.set_use_toroidal_cores(false);
        settings.set_use_concentric_cores(true);
        json inputsJson;


        auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "psim_simulation.csv").string();

        double frequency = 120000;
        auto reader = CircuitSimulationReader(simulation_path);
        auto operatingPoint = reader.extract_operating_point(2, frequency);

        operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 52e-6);


        inputsJson["designRequirements"] = json();
        inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 52e-6;
        inputsJson["designRequirements"]["turnsRatios"] = json::array();
        json turnRatioSecondary;
        turnRatioSecondary["minimum"] = 2;
        turnRatioSecondary["maximum"] = 2.5;
        turnRatioSecondary["nominal"] = 2.275;
        inputsJson["designRequirements"]["turnsRatios"].push_back(turnRatioSecondary);
        inputsJson["operatingPoints"] = json::array();
        json operatingPointJson;
        to_json(operatingPointJson, operatingPoint);
        inputsJson["operatingPoints"].push_back(operatingPointJson);

        OpenMagnetics::Inputs inputs(inputsJson);


        std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY,
                                                                    IsolationSide::SECONDARY};


        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process();
        MagneticAdviser MagneticAdviser;

        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 3);

        REQUIRE(masMagnetics.size() > 0);

        if (plot)
        for (auto [masMagnetic, scoring] : masMagnetics) {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");

            {
                auto outFile = outputFilePath;
                std::string filename = "Test_Magnetic_Adviser_Psim_" +  std::to_string(scoring) + ".mas.json";
                outFile.append(filename);
                to_file(outFile, masMagnetic);
            }

            // MagneticAdviser::preview_magnetic(masMagnetic);
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outFile = outputFilePath;
            std::string filename = "Test_Magnetic_Adviser_Psim_" + std::to_string(scoring) + ".svg";
            outFile.append(filename);
            Painter painter(outFile, true);

            settings.set_painter_mode(PainterModes::CONTOUR);
            settings.set_painter_logarithmic_scale(true);
            settings.set_painter_include_fringing(true);
            settings.set_painter_maximum_value_colorbar(std::nullopt);
            settings.set_painter_minimum_value_colorbar(std::nullopt);
            // settings.set_painter_number_points_x(200);
            // settings.set_painter_number_points_y(200);
            painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
            painter.paint_core(masMagnetic.get_mutable_magnetic());
            // painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            #ifdef ENABLE_MATPLOTPP
            painter.export_svg();
            #else
                INFO("matplotplusplus disabled — skipping AdvancedPainter call");
            #endif

        }
    }

    TEST_CASE("Test_Magnetic_Adviser_Plecs", "[adviser][magnetic-adviser]") {
        settings.reset();
        clear_databases();
        settings.set_use_only_cores_in_stock(true);
        settings.set_use_toroidal_cores(false);
        settings.set_use_concentric_cores(true);
        json inputsJson;


        auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "plecs_simulation.csv").string();

        double frequency = 50;
        auto reader = CircuitSimulationReader(simulation_path);
        auto operatingPoint = reader.extract_operating_point(1, frequency);

        operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 100e-6);


        inputsJson["designRequirements"] = json();
        inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
        inputsJson["designRequirements"]["turnsRatios"] = json::array();
        inputsJson["operatingPoints"] = json::array();
        json operatingPointJson;
        to_json(operatingPointJson, operatingPoint);
        inputsJson["operatingPoints"].push_back(operatingPointJson);

        OpenMagnetics::Inputs inputs(inputsJson);


        std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY};


        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process();
        MagneticAdviser MagneticAdviser;

        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 3);

        REQUIRE(masMagnetics.size() > 0);

        if (plot)
        for (auto [masMagnetic, scoring] : masMagnetics) {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");

            {
                auto outFile = outputFilePath;
                std::string filename = "Test_Magnetic_Adviser_Plecs_" +  std::to_string(scoring) + ".mas.json";
                outFile.append(filename);
                to_file(outFile, masMagnetic);
            }

            // MagneticAdviser::preview_magnetic(masMagnetic);
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outFile = outputFilePath;
            std::string filename = "Test_Magnetic_Adviser_Plecs_" + std::to_string(scoring) + ".svg";
            outFile.append(filename);
            Painter painter(outFile, true);

            settings.set_painter_mode(PainterModes::CONTOUR);
            settings.set_painter_logarithmic_scale(true);
            settings.set_painter_include_fringing(true);
            settings.set_painter_maximum_value_colorbar(std::nullopt);
            settings.set_painter_minimum_value_colorbar(std::nullopt);
            // settings.set_painter_number_points_x(200);
            // settings.set_painter_number_points_y(200);
            painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
            painter.paint_core(masMagnetic.get_mutable_magnetic());
            // painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            #ifdef ENABLE_MATPLOTPP
            painter.export_svg();
            #else
                INFO("matplotplusplus disabled — skipping AdvancedPainter call");
            #endif

        }
    }

    TEST_CASE("Test_Magnetic_Adviser_Ltspice", "[adviser][magnetic-adviser]") {
        SKIP("Test needs investigation");
        settings.reset();
        clear_databases();
        settings.set_use_only_cores_in_stock(true);
        settings.set_use_toroidal_cores(false);
        settings.set_use_concentric_cores(true);
        json inputsJson;


        auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "ltspice_simulation.txt").string();

        double frequency = 372618;
        auto reader = CircuitSimulationReader(simulation_path);
        auto operatingPoint = reader.extract_operating_point(2, frequency);

        operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 100e-6);


        inputsJson["designRequirements"] = json();
        inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
        inputsJson["designRequirements"]["turnsRatios"] = json::array();
        json turnRatioSecondary;
        turnRatioSecondary["nominal"] = 1.0 / 7.11;
        inputsJson["designRequirements"]["turnsRatios"].push_back(turnRatioSecondary);
        inputsJson["operatingPoints"] = json::array();
        json operatingPointJson;
        to_json(operatingPointJson, operatingPoint);
        inputsJson["operatingPoints"].push_back(operatingPointJson);

        OpenMagnetics::Inputs inputs(inputsJson);


        std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY,
                                                                    IsolationSide::SECONDARY};


        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process();
        MagneticAdviser MagneticAdviser;

        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 3);

        REQUIRE(masMagnetics.size() > 0);

        if (plot)
        for (auto [masMagnetic, scoring] : masMagnetics) {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");

            {
                auto outFile = outputFilePath;
                std::string filename = "Test_Magnetic_Adviser_Ltspice_" +  std::to_string(scoring) + ".mas.json";
                outFile.append(filename);
                to_file(outFile, masMagnetic);
            }

            // MagneticAdviser::preview_magnetic(masMagnetic);
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outFile = outputFilePath;
            std::string filename = "Test_Magnetic_Adviser_Ltspice_" + std::to_string(scoring) + ".svg";
            outFile.append(filename);
            Painter painter(outFile, true);

            settings.set_painter_mode(PainterModes::CONTOUR);
            settings.set_painter_logarithmic_scale(true);
            settings.set_painter_include_fringing(true);
            settings.set_painter_maximum_value_colorbar(std::nullopt);
            settings.set_painter_minimum_value_colorbar(std::nullopt);
            // settings.set_painter_number_points_x(200);
            // settings.set_painter_number_points_y(200);
            painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
            painter.paint_core(masMagnetic.get_mutable_magnetic());
            // painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            #ifdef ENABLE_MATPLOTPP
            painter.export_svg();
            #else
                INFO("matplotplusplus disabled — skipping AdvancedPainter call");
            #endif

        }
    }

    TEST_CASE("Test_Magnetic_Adviser_Simba", "[adviser][magnetic-adviser]") {
        settings.reset();
        clear_databases();
        settings.set_use_only_cores_in_stock(true);
        settings.set_use_toroidal_cores(false);
        settings.set_use_concentric_cores(true);
        json inputsJson;


        auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "simba_simulation.csv").string();

        double frequency = 100000;
        auto reader = CircuitSimulationReader(simulation_path);
        auto operatingPoint = reader.extract_operating_point(2, frequency);

        operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 220e-6);


        inputsJson["designRequirements"] = json();
        inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 220e-6;
        inputsJson["designRequirements"]["turnsRatios"] = json::array();
        json turnRatioSecondary;
        turnRatioSecondary["minimum"] = 3;
        turnRatioSecondary["nominal"] = 3.33;
        turnRatioSecondary["maximum"] = 3.5;
        inputsJson["designRequirements"]["turnsRatios"].push_back(turnRatioSecondary);
        inputsJson["operatingPoints"] = json::array();
        json operatingPointJson;

        to_json(operatingPointJson, operatingPoint);
        inputsJson["operatingPoints"].push_back(operatingPointJson);

        OpenMagnetics::Inputs inputs(inputsJson);


        std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY,
                                                                    IsolationSide::SECONDARY};


        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process();
        MagneticAdviser MagneticAdviser;

        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 3);

        REQUIRE(masMagnetics.size() > 0);

        if (plot)
        for (auto [masMagnetic, scoring] : masMagnetics) {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");

            {
                auto outFile = outputFilePath;
                std::string filename = "Test_Magnetic_Adviser_Ltspice_" +  std::to_string(scoring) + ".mas.json";
                outFile.append(filename);
                to_file(outFile, masMagnetic);
            }

            // MagneticAdviser::preview_magnetic(masMagnetic);
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outFile = outputFilePath;
            std::string filename = "Test_Magnetic_Adviser_Ltspice_" + std::to_string(scoring) + ".svg";
            outFile.append(filename);
            Painter painter(outFile, true);

            settings.set_painter_mode(PainterModes::CONTOUR);
            settings.set_painter_logarithmic_scale(true);
            settings.set_painter_include_fringing(true);
            settings.set_painter_maximum_value_colorbar(std::nullopt);
            settings.set_painter_minimum_value_colorbar(std::nullopt);
            // settings.set_painter_number_points_x(200);
            // settings.set_painter_number_points_y(200);
            painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
            painter.paint_core(masMagnetic.get_mutable_magnetic());
            // painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            #ifdef ENABLE_MATPLOTPP
            painter.export_svg();
            #else
                INFO("matplotplusplus disabled — skipping AdvancedPainter call");
            #endif

        }
    }

    TEST_CASE("Test_Magnetic_Adviser_Rosano_Forward", "[adviser][magnetic-adviser]") {
        settings.reset();
        clear_databases();
        settings.set_use_only_cores_in_stock(false);
        settings.set_use_toroidal_cores(false);
        settings.set_use_concentric_cores(true);
        json inputsJson;


        auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "forward_case.csv").string();

        double frequency = 200000;
        auto reader = CircuitSimulationReader(simulation_path);
        std::vector<std::map<std::string, std::string>> mapColumnNames;
        std::map<std::string, std::string> primaryColumnNames;
        std::map<std::string, std::string> secondaryColumnNames;
        primaryColumnNames["time"] = "Time";
        primaryColumnNames["current"] = "Ipri";
        primaryColumnNames["voltage"] = "Vpri";
        mapColumnNames.push_back(primaryColumnNames);
        secondaryColumnNames["time"] = "Time";
        secondaryColumnNames["current"] = "Isec";
        secondaryColumnNames["voltage"] = "Vout";
        mapColumnNames.push_back(secondaryColumnNames);
        auto operatingPoint = reader.extract_operating_point(2, frequency, mapColumnNames);

        operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 121e-6, std::vector<double>({1.33}));


        inputsJson["designRequirements"] = json();
        inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 121e-6;
        inputsJson["designRequirements"]["turnsRatios"] = json::array();
        json turnRatioSecondary;
        turnRatioSecondary["minimum"] = 1.3;
        turnRatioSecondary["nominal"] = 1.33;
        turnRatioSecondary["maximum"] = 1.4;
        inputsJson["designRequirements"]["turnsRatios"].push_back(turnRatioSecondary);
        inputsJson["operatingPoints"] = json::array();
        json operatingPointJson;

        to_json(operatingPointJson, operatingPoint);
        inputsJson["operatingPoints"].push_back(operatingPointJson);

        OpenMagnetics::Inputs inputs(inputsJson);


        std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY,
                                                                    IsolationSide::SECONDARY};


        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process();
        MagneticAdviser MagneticAdviser;

        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 3);

        REQUIRE(masMagnetics.size() > 0);

        if (plot)
        for (auto [masMagnetic, scoring] : masMagnetics) {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");

            {
                auto outFile = outputFilePath;
                std::string filename = "Test_Magnetic_Adviser_Rosano_Forward_" +  std::to_string(scoring) + ".mas.json";
                outFile.append(filename);
                to_file(outFile, masMagnetic);
            }

            // MagneticAdviser::preview_magnetic(masMagnetic);
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outFile = outputFilePath;
            std::string filename = "Test_Magnetic_Adviser_Rosano_Forward_" + std::to_string(scoring) + ".svg";
            outFile.append(filename);
            Painter painter(outFile, true);

            settings.set_painter_mode(PainterModes::CONTOUR);
            settings.set_painter_logarithmic_scale(true);
            settings.set_painter_include_fringing(true);
            settings.set_painter_maximum_value_colorbar(std::nullopt);
            settings.set_painter_minimum_value_colorbar(std::nullopt);
            // settings.set_painter_number_points_x(200);
            // settings.set_painter_number_points_y(200);
            painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
            painter.paint_core(masMagnetic.get_mutable_magnetic());
            // painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            #ifdef ENABLE_MATPLOTPP
            painter.export_svg();
            #else
                INFO("matplotplusplus disabled — skipping AdvancedPainter call");
            #endif

        }
    }

    TEST_CASE("Test_Magnetic_Adviser_Rosano_Flyback", "[adviser][magnetic-adviser]") {
        settings.reset();
        clear_databases();
        settings.set_use_only_cores_in_stock(false);
        settings.set_use_toroidal_cores(false);
        settings.set_use_concentric_cores(true);
        json inputsJson;


        auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "flyback_case.csv").string();

        double frequency = 200000;
        auto reader = CircuitSimulationReader(simulation_path);
        std::vector<std::map<std::string, std::string>> mapColumnNames;
        std::map<std::string, std::string> primaryColumnNames;
        std::map<std::string, std::string> secondaryColumnNames;
        primaryColumnNames["time"] = "Time";
        primaryColumnNames["current"] = "Ipri";
        primaryColumnNames["voltage"] = "Vpri";
        mapColumnNames.push_back(primaryColumnNames);
        secondaryColumnNames["time"] = "Time";
        secondaryColumnNames["current"] = "Isec";
        secondaryColumnNames["voltage"] = "Vsec";
        mapColumnNames.push_back(secondaryColumnNames);
        auto operatingPoint = reader.extract_operating_point(2, frequency, mapColumnNames);

        operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 121e-6, std::vector<double>({1.33}));


        inputsJson["designRequirements"] = json();
        inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 121e-6;
        inputsJson["designRequirements"]["turnsRatios"] = json::array();
        json turnRatioSecondary;
        turnRatioSecondary["minimum"] = 1.3;
        turnRatioSecondary["nominal"] = 1.33;
        turnRatioSecondary["maximum"] = 1.4;
        inputsJson["designRequirements"]["turnsRatios"].push_back(turnRatioSecondary);
        inputsJson["operatingPoints"] = json::array();
        json operatingPointJson;

        to_json(operatingPointJson, operatingPoint);
        inputsJson["operatingPoints"].push_back(operatingPointJson);

        OpenMagnetics::Inputs inputs(inputsJson);


        std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY,
                                                                    IsolationSide::SECONDARY};


        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process();
        MagneticAdviser MagneticAdviser;

        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 3);

        REQUIRE(masMagnetics.size() > 0);

        if (plot)
        for (auto [masMagnetic, scoring] : masMagnetics) {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");

            {
                auto outFile = outputFilePath;
                std::string filename = "Test_Magnetic_Adviser_Rosano_Forward_" +  std::to_string(scoring) + ".mas.json";
                outFile.append(filename);
                to_file(outFile, masMagnetic);
            }

            // MagneticAdviser::preview_magnetic(masMagnetic);
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outFile = outputFilePath;
            std::string filename = "Test_Magnetic_Adviser_Rosano_Forward_" + std::to_string(scoring) + ".svg";
            outFile.append(filename);
            Painter painter(outFile, true);

            settings.set_painter_mode(PainterModes::CONTOUR);
            settings.set_painter_logarithmic_scale(true);
            settings.set_painter_include_fringing(true);
            settings.set_painter_maximum_value_colorbar(std::nullopt);
            settings.set_painter_minimum_value_colorbar(std::nullopt);
            // settings.set_painter_number_points_x(200);
            // settings.set_painter_number_points_y(200);
            painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
            painter.paint_core(masMagnetic.get_mutable_magnetic());
            // painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            #ifdef ENABLE_MATPLOTPP
            painter.export_svg();
            #else
                INFO("matplotplusplus disabled — skipping AdvancedPainter call");
            #endif

        }
    }

    TEST_CASE("Test_Magnetic_Adviser_PFC_Only_Current", "[adviser][magnetic-adviser]") {
        settings.reset();
        clear_databases();
        settings.set_use_only_cores_in_stock(true);
        settings.set_use_toroidal_cores(false);
        settings.set_use_concentric_cores(true);
        json inputsJson;


        auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "only_pfc_current_waveform.csv").string();

        double frequency = 50;
        auto reader = CircuitSimulationReader(simulation_path);
        auto operatingPoint = reader.extract_operating_point(1, frequency);

        operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 110e-6);


        inputsJson["designRequirements"] = json();
        inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 110e-6;
        inputsJson["designRequirements"]["turnsRatios"] = json::array();
        inputsJson["operatingPoints"] = json::array();
        json operatingPointJson;

        to_json(operatingPointJson, operatingPoint);
        inputsJson["operatingPoints"].push_back(operatingPointJson);

        OpenMagnetics::Inputs inputs(inputsJson);


        std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY};


        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process();
        MagneticAdviser MagneticAdviser;

        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 3);

        REQUIRE(masMagnetics.size() > 0);

        if (plot)
        for (auto [masMagnetic, scoring] : masMagnetics) {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");

            {
                auto outFile = outputFilePath;
                std::string filename = "Test_Magnetic_Adviser_PFC_Only_Current_" +  std::to_string(scoring) + ".mas.json";
                outFile.append(filename);
                to_file(outFile, masMagnetic);
            }

            // MagneticAdviser::preview_magnetic(masMagnetic);
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outFile = outputFilePath;
            std::string filename = "Test_Magnetic_Adviser_PFC_Only_Current_" + std::to_string(scoring) + ".svg";
            outFile.append(filename);
            Painter painter(outFile, true);

            settings.set_painter_mode(PainterModes::CONTOUR);
            settings.set_painter_logarithmic_scale(true);
            settings.set_painter_include_fringing(true);
            settings.set_painter_maximum_value_colorbar(std::nullopt);
            settings.set_painter_minimum_value_colorbar(std::nullopt);
            // settings.set_painter_number_points_x(200);
            // settings.set_painter_number_points_y(200);
            painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
            painter.paint_core(masMagnetic.get_mutable_magnetic());
            // painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            #ifdef ENABLE_MATPLOTPP
            painter.export_svg();
            #else
                INFO("matplotplusplus disabled — skipping AdvancedPainter call");
            #endif

        }
    }

    TEST_CASE("Test_Magnetic_Adviser_CMC", "[adviser][magnetic-adviser]") {
        clear_databases();
        settings.set_use_concentric_cores(false);
        settings.set_use_toroidal_cores(true);
        settings.set_use_only_cores_in_stock(false);
        double voltagePeakToPeak = 600;
        double dcCurrent = 30;
        double ambientTemperature = 25;
        double frequency = 100000;
        double desiredMagnetizingInductance = 10e-5;
        std::vector<double> turnsRatios = {};
        OpenMagnetics::Inputs inputs;


        DesignRequirements designRequirements = inputs.get_design_requirements();
        inputs = OpenMagnetics::Inputs::create_quick_operating_point(frequency, desiredMagnetizingInductance, ambientTemperature, WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, dcCurrent, turnsRatios);
        if (designRequirements.get_insulation()) {
            inputs.get_mutable_design_requirements().set_insulation(designRequirements.get_insulation().value());
        }

        std::vector<std::pair<double, double>> impedancePoints = {
            {1e6, 1000},
            {2e6, 2000},
        };

        std::vector<ImpedanceAtFrequency> minimumImpedance;
        for (auto [frequencyPoint, impedanceMagnitudePoint] : impedancePoints) {
            ImpedancePoint impedancePoint;
            impedancePoint.set_magnitude(impedanceMagnitudePoint);
            ImpedanceAtFrequency impedanceAtFrequency;
            impedanceAtFrequency.set_frequency(frequencyPoint);
            impedanceAtFrequency.set_impedance(impedancePoint);
            minimumImpedance.push_back(impedanceAtFrequency);
        }
        inputs.get_mutable_design_requirements().set_minimum_impedance(minimumImpedance);

        MagneticAdviser magneticAdviser;
        auto masMagnetics = magneticAdviser.get_advised_magnetic(inputs, 50);


        // REQUIRE(masMagnetics.size() == 1);
        // double bestScoring = masMagnetics[0].second;
        // for (size_t i = 0; i < masMagnetics.size(); ++i)
        // {
        //     REQUIRE(masMagnetics[i].second <= bestScoring);
        // }

        // REQUIRE(masMagnetics[0].first.get_magnetic().get_core().get_name() == "T 25.3/14.8/20 - N30 - Ungapped");
        // auto magnetic = masMagnetics[0].first.get_magnetic();

        // auto selfResonantFrequencyFast = Impedance().calculate_self_resonant_frequency(masMagnetics[0].first.get_magnetic());

        // for (auto [frequencyPoint, impedanceMagnitudePoint] : impedancePoints) {
        //     auto impedance = Impedance().calculate_impedance(masMagnetics[0].first.get_magnetic(), frequencyPoint);
        //     REQUIRE(frequencyPoint < selfResonantFrequencyFast * 0.50);
        //     REQUIRE(abs(impedance) >= impedanceMagnitudePoint);
        // }

        // auto selfResonantFrequency = Impedance(false).calculate_self_resonant_frequency(masMagnetics[0].first.get_magnetic());
        // for (auto [frequencyPoint, impedanceMagnitudePoint] : impedancePoints) {
        //     auto impedance = Impedance(false).calculate_impedance(masMagnetics[0].first.get_magnetic(), frequencyPoint);
        //     REQUIRE(frequencyPoint < selfResonantFrequency * 0.50);
        //     REQUIRE(abs(impedance) >= impedanceMagnitudePoint);
        // }

        // {
        //     auto impedanceSweep = Sweeper().sweep_impedance_over_frequency(masMagnetics[0].first.get_magnetic(), 1000, 40000000, 1000);

        //     auto outputFilePath = std::filesystem::path {std::source_location::current().file_name()}.parent_path().append("..").append("output");
        //     auto outFile = outputFilePath;

        //     outFile.append("Test_Magnetic_Adviser_CMC_Impedance.svg");
        //     std::filesystem::remove(outFile);
        //     Painter painter(outFile, false, true);
        //     painter.paint_curve(impedanceSweep, true);
        //     painter.export_svg();
        //     REQUIRE(std::filesystem::exists(outFile));
        // }
        {
            auto outputFilePath = std::filesystem::path {std::source_location::current().file_name()}.parent_path().append("..").append("output");
            auto outFile = outputFilePath;

            outFile.append("Test_Magnetic_Adviser_CMC.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_core(masMagnetics[0].first.get_magnetic());
            painter.paint_bobbin(masMagnetics[0].first.get_magnetic());
            painter.paint_coil_turns(masMagnetics[0].first.get_magnetic());
            #ifdef ENABLE_MATPLOTPP
            painter.export_svg();
            #else
                INFO("matplotplusplus disabled — skipping AdvancedPainter call");
            #endif

            REQUIRE(std::filesystem::exists(outFile));

        }

        settings.reset();
    }

    TEST_CASE("Test_CatalogueAdviser_Found", "[adviser][magnetic-adviser][catalogue]") { 

        {
            auto external_core_materials_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "core_materials.ndjson").string();

            std::ifstream file(external_core_materials_path, std::ios_base::binary | std::ios_base::in);
            if(!file.is_open())
                throw std::runtime_error("Failed to open " + external_core_materials_path);
            using Iterator = std::istreambuf_iterator<char>;
            std::string external_core_materials(Iterator{file}, Iterator{});

            clear_databases();
            load_core_materials();

            auto allCoreMaterials = coreMaterialDatabase;

            load_core_materials(external_core_materials);
        }

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {24, 24};
        std::vector<int64_t> numberParallels = {1, 1};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 107026;
        double magnetizingInductance = 100e-6;
        double temperature = 25;
        WaveformLabel waveShape = WaveformLabel::TRIANGULAR;
        double peakToPeak = 20;
        double dutyCycle = 0.5;
        double dcCurrent = 0;

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);


        auto requirements = inputs.get_mutable_design_requirements();
        auto insulationRequirements = requirements.get_insulation().value();
        auto standards = std::vector<InsulationStandards>{InsulationStandards::IEC_606641};
        insulationRequirements.set_standards(standards);
        requirements.set_insulation(insulationRequirements);
        MaximumDimensions maximumDimensions;
        maximumDimensions.set_depth(0.05);
        maximumDimensions.set_height(0.05);
        maximumDimensions.set_width(0.05);
        requirements.set_maximum_dimensions(maximumDimensions);

        DimensionWithTolerance magnetizingInductanceWithTolerance;
        magnetizingInductanceWithTolerance.set_minimum(magnetizingInductance);
        magnetizingInductanceWithTolerance.set_nominal(std::nullopt);
        magnetizingInductanceWithTolerance.set_maximum(std::nullopt);

        requirements.set_magnetizing_inductance(magnetizingInductanceWithTolerance);

        inputs.set_design_requirements(requirements);

        OpenMagnetics::Mas masMagnetic;
        inputs.process();


        std::vector<std::pair<double, double>> impedancePoints = {
            {1e5, 10000},
            {2e5, 20000},
        };

        std::vector<ImpedanceAtFrequency> minimumImpedance;
        for (auto [frequencyPoint, impedanceMagnitudePoint] : impedancePoints) {
            ImpedancePoint impedancePoint;
            impedancePoint.set_magnitude(impedanceMagnitudePoint);
            ImpedanceAtFrequency impedanceAtFrequency;
            impedanceAtFrequency.set_frequency(frequencyPoint);
            impedanceAtFrequency.set_impedance(impedancePoint);
            minimumImpedance.push_back(impedanceAtFrequency);
        }
        inputs.get_mutable_design_requirements().set_minimum_impedance(minimumImpedance);



        std::vector<OpenMagnetics::Magnetic> catalogue;
        auto inventory_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "cmcs.ndjson").string();

        std::ifstream ndjsonFile(inventory_path);
        std::string jsonLine;
        json coresJson = json::array();
        while (std::getline(ndjsonFile, jsonLine)) {
            json jf = json::parse(jsonLine);
            OpenMagnetics::Magnetic magnetic(jf);
            catalogue.push_back(magnetic);
        }


        MagneticAdviser magneticAdviser;
        auto masMagnetics = magneticAdviser.get_advised_magnetic(inputs, catalogue, 1);
        REQUIRE(masMagnetics.size() > 0);

        auto scorings = magneticAdviser.get_scorings();
        for (auto [name, values] : scorings) {
            double scoringTotal = 0;
            for (auto [key, value] : values) {
                scoringTotal += value;
            }
        }

        for (auto masMagneticWithScoring : masMagnetics) {
            auto masMagnetic = masMagneticWithScoring.first;


            auto operatingPoint = inputs.get_operating_points()[0];
            MagnetizingInductance magnetizingInductanceObj;
            auto magneticFluxDensity = magnetizingInductanceObj.calculate_inductance_and_magnetic_flux_density(masMagnetic.get_mutable_magnetic().get_core(), masMagnetic.get_mutable_magnetic().get_coil(), &operatingPoint).second;

            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            if (plot) {
                auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
                auto outFile = outputFilePath;
                std::string filename = "MagneticAdviser" + masMagnetic.get_magnetic().get_manufacturer_info().value().get_reference().value() + ".svg";
                outFile.append(filename);
                Painter painter(outFile, true);

                painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
                painter.paint_core(masMagnetic.get_mutable_magnetic());
                painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
                painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
                #ifdef ENABLE_MATPLOTPP
                painter.export_svg();
                #else
                    INFO("matplotplusplus disabled — skipping AdvancedPainter call");
                #endif

            }
        }
    }

    TEST_CASE("Test_CatalogueAdviser_Not_Found", "[adviser][magnetic-adviser][catalogue]") {
        {
            auto external_core_materials_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "core_materials.ndjson").string();

            std::ifstream file(external_core_materials_path, std::ios_base::binary | std::ios_base::in);
            if(!file.is_open())
                throw std::runtime_error("Failed to open " + external_core_materials_path);
            using Iterator = std::istreambuf_iterator<char>;
            std::string external_core_materials(Iterator{file}, Iterator{});

            clear_databases();
            load_core_materials();

            auto allCoreMaterials = coreMaterialDatabase;

            load_core_materials(external_core_materials);
        }

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {24, 24};
        std::vector<int64_t> numberParallels = {1, 1};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 507026;
        double magnetizingInductance = 100e-6;
        double temperature = 25;
        WaveformLabel waveShape = WaveformLabel::TRIANGULAR;
        double peakToPeak = 1;
        double dutyCycle = 0.5;
        double dcCurrent = 0;

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);


        auto requirements = inputs.get_mutable_design_requirements();
        auto insulationRequirements = requirements.get_insulation().value();
        auto standards = std::vector<InsulationStandards>{InsulationStandards::IEC_606641};
        insulationRequirements.set_standards(standards);
        requirements.set_insulation(insulationRequirements);
        MaximumDimensions maximumDimensions;
        maximumDimensions.set_depth(0.001);
        maximumDimensions.set_height(0.001);
        maximumDimensions.set_width(0.001);
        requirements.set_maximum_dimensions(maximumDimensions);

        DimensionWithTolerance magnetizingInductanceWithTolerance;
        magnetizingInductanceWithTolerance.set_minimum(magnetizingInductance);
        magnetizingInductanceWithTolerance.set_nominal(std::nullopt);
        magnetizingInductanceWithTolerance.set_maximum(std::nullopt);

        requirements.set_magnetizing_inductance(magnetizingInductanceWithTolerance);

        inputs.set_design_requirements(requirements);

        OpenMagnetics::Mas masMagnetic;
        inputs.process();


        std::vector<std::pair<double, double>> impedancePoints = {
            {1e5, 10000},
            {2e5, 20000},
        };

        std::vector<ImpedanceAtFrequency> minimumImpedance;
        for (auto [frequencyPoint, impedanceMagnitudePoint] : impedancePoints) {
            ImpedancePoint impedancePoint;
            impedancePoint.set_magnitude(impedanceMagnitudePoint);
            ImpedanceAtFrequency impedanceAtFrequency;
            impedanceAtFrequency.set_frequency(frequencyPoint);
            impedanceAtFrequency.set_impedance(impedancePoint);
            minimumImpedance.push_back(impedanceAtFrequency);
        }
        inputs.get_mutable_design_requirements().set_minimum_impedance(minimumImpedance);



        std::vector<OpenMagnetics::Magnetic> catalogue;
        auto inventory_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "cmcs.ndjson").string();

        std::ifstream ndjsonFile(inventory_path);
        std::string jsonLine;
        json coresJson = json::array();
        while (std::getline(ndjsonFile, jsonLine)) {
            json jf = json::parse(jsonLine);
            OpenMagnetics::Magnetic magnetic(jf);
            catalogue.push_back(magnetic);
        }


        MagneticAdviser magneticAdviser;
        auto masMagnetics = magneticAdviser.get_advised_magnetic(inputs, catalogue, 1);
        REQUIRE(masMagnetics.size() > 0);

        auto scorings = magneticAdviser.get_scorings();
        for (auto [name, values] : scorings) {
            double scoringTotal = 0;
            for (auto [key, value] : values) {
                scoringTotal += value;
            }
        }

        for (auto masMagneticWithScoring : masMagnetics) {
            auto masMagnetic = masMagneticWithScoring.first;
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            if (plot) {
                auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
                auto outFile = outputFilePath;
                std::string filename = "MagneticAdviser" + masMagnetic.get_magnetic().get_manufacturer_info().value().get_reference().value() + ".svg";
                outFile.append(filename);
                Painter painter(outFile, true);

                painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
                painter.paint_core(masMagnetic.get_mutable_magnetic());
                painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
                painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
                #ifdef ENABLE_MATPLOTPP
                painter.export_svg();
                #else
                    INFO("matplotplusplus disabled — skipping AdvancedPainter call");
                #endif

            }
        }
    }

    TEST_CASE("Test_CatalogueAdviser_Web_0", "[adviser][magnetic-adviser][catalogue][bug]") {
        SKIP("Test needs investigation");

        {
            auto external_core_materials_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "core_materials.ndjson").string();

            std::ifstream file(external_core_materials_path, std::ios_base::binary | std::ios_base::in);
            if(!file.is_open())
                throw std::runtime_error("Failed to open " + external_core_materials_path);
            using Iterator = std::istreambuf_iterator<char>;
            std::string external_core_materials(Iterator{file}, Iterator{});

            clear_databases();
            load_core_materials();

            auto allCoreMaterials = coreMaterialDatabase;

            load_core_materials(external_core_materials);
        }

        OpenMagnetics::Inputs inputs;
        {
            auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "bug_catalogue.json").string();
            auto mas = OpenMagneticsTesting::mas_loader(path);
            inputs = mas.get_inputs();
        }

        std::vector<OpenMagnetics::Magnetic> catalogue;
        {
            auto inventory_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "cmcs.ndjson").string();

            std::ifstream ndjsonFile(inventory_path);
            std::string jsonLine;
            json coresJson = json::array();
            while (std::getline(ndjsonFile, jsonLine)) {
                json jf = json::parse(jsonLine);
                OpenMagnetics::Magnetic magnetic(jf);
                catalogue.push_back(magnetic);
            }
        }


        MagneticAdviser magneticAdviser;
        auto masMagnetics = magneticAdviser.get_advised_magnetic(inputs, catalogue, 1);
        REQUIRE(masMagnetics.size() > 0);

        auto scorings = magneticAdviser.get_scorings();
        for (auto [name, values] : scorings) {
            double scoringTotal = 0;
            for (auto [key, value] : values) {
                scoringTotal += value;
            }
        }

        for (auto masMagneticWithScoring : masMagnetics) {
            auto masMagnetic = masMagneticWithScoring.first;
            auto magnetic = masMagnetic.get_mutable_magnetic();
            auto core = magnetic.get_mutable_core();
            auto coil = magnetic.get_mutable_coil();
            auto operatingPoint = masMagnetic.get_inputs().get_operating_points()[0];

            for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
                auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
                if (!excitation.get_current()) {
                    throw std::runtime_error("Current is missing in excitation");
                }
            }
            OpenMagneticsTesting::check_turns_description(magnetic.get_coil());
            if (plot) {
                auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
                auto outFile = outputFilePath;
                std::string filename = "MagneticAdviser" + masMagnetic.get_magnetic().get_manufacturer_info().value().get_reference().value() + ".svg";
                outFile.append(filename);
                Painter painter(outFile, true);

                painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
                painter.paint_core(masMagnetic.get_mutable_magnetic());
                painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
                painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
                #ifdef ENABLE_MATPLOTPP
                painter.export_svg();
                #else
                    INFO("matplotplusplus disabled — skipping AdvancedPainter call");
                #endif

            }
        }
    }

    TEST_CASE("Test_CatalogueAdviser_Web_1", "[adviser][magnetic-adviser][catalogue][bug][smoke-test]") {

        auto json_path_2218 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_catalogueadviser_web_1_2218.json");
        std::ifstream json_file_2218(json_path_2218);
        json catalogueJson = json::parse(json_file_2218);

        std::vector<OpenMagnetics::Mas> catalogue;
        from_json(catalogueJson, catalogue);

        json filterFlowJson = json::parse(R"([{"filter": "Dc Current Density", "invert": false, "log": false, "strictlyRequired": false, "weight": 1 } ])");
        std::vector<OpenMagnetics::MagneticFilterOperation> filterFlow = filterFlowJson.get<std::vector<OpenMagnetics::MagneticFilterOperation>>();

        MagneticAdviser magneticAdviser(false);
        auto masMagnetics = magneticAdviser.get_advised_magnetic(catalogue, filterFlow, 2);
        REQUIRE(masMagnetics.size() > 0);

        auto scorings = magneticAdviser.get_scorings();
        for (auto [name, values] : scorings) {
            double scoringTotal = 0;
            for (auto [key, value] : values) {
                scoringTotal += value;
            }
        }
    }

    TEST_CASE("Test_CatalogueAdviser_Web_2", "[adviser][magnetic-adviser][catalogue][bug][smoke-test]") {

        auto json_path_2241 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_catalogueadviser_web_2_2241.json");
        std::ifstream json_file_2241(json_path_2241);
        json catalogueJson = json::parse(json_file_2241);

        std::vector<OpenMagnetics::Mas> catalogue;
        from_json(catalogueJson, catalogue);

        json filterFlowJson = json::parse(R"([{"filter": "Turns Ratios", "invert": false, "log": false, "strictlyRequired": true, "weight": 0.1 }, {"filter": "Magnetizing Inductance", "invert": false, "log": false, "strictlyRequired": false, "weight": 0.5 }, {"filter": "Core And DC Losses", "invert": false, "log": false, "strictlyRequired": false, "weight": 1 } ])");
        std::vector<OpenMagnetics::MagneticFilterOperation> filterFlow = filterFlowJson.get<std::vector<OpenMagnetics::MagneticFilterOperation>>();

        MagneticAdviser magneticAdviser(false);
        auto masMagnetics = magneticAdviser.get_advised_magnetic(catalogue, filterFlow, 2);
        REQUIRE(masMagnetics.size() > 0);
    }

    TEST_CASE("Test_CatalogueAdviser_Web_3", "[adviser][magnetic-adviser][catalogue][bug]") {

        auto json_path_2256 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_catalogueadviser_web_3_2256.json");
        std::ifstream json_file_2256(json_path_2256);
        json catalogueJson = json::parse(json_file_2256);

        std::vector<OpenMagnetics::Mas> catalogue;
        from_json(catalogueJson, catalogue);

        json filterFlowJson = json::parse(R"([{"filter": "Turns Ratios", "invert": false, "log": false, "strictlyRequired": true, "weight": 0.1 }, {"filter": "Magnetizing Inductance", "invert": false, "log": false, "strictlyRequired": false, "weight": 0.5 }, {"filter": "Losses", "invert": false, "log": false, "strictlyRequired": false, "weight": 1 } ])");
        std::vector<OpenMagnetics::MagneticFilterOperation> filterFlow = filterFlowJson.get<std::vector<OpenMagnetics::MagneticFilterOperation>>();

        MagneticAdviser magneticAdviser(false);
        auto masMagnetics = magneticAdviser.get_advised_magnetic(catalogue, filterFlow, 2);
        REQUIRE(masMagnetics.size() > 0);
    }

    TEST_CASE("Test_CatalogueAdviser_Web_4", "[adviser][magnetic-adviser][catalogue][bug]") {
        SKIP("Test needs investigation");

        auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "catalogue_missing_turns.json").string();
        auto mas = OpenMagneticsTesting::mas_loader(path);

        auto catalogue_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "magnetics.ndjson").string();

        std::ifstream ndjsonFile(catalogue_path);
        std::string jsonLine;
        json coresJson = json::array();
        while (std::getline(ndjsonFile, jsonLine)) {
            json jf = json::parse(jsonLine);
            OpenMagnetics::Magnetic magnetic(jf); 
            magneticsCache.load(magnetic.get_reference(), magnetic);
        }
        magneticsCache.autocomplete_magnetics();

        json filterFlowJson = json::parse(R"([{"filter": "Turns Ratios", "invert": true, "log": false, "strictlyRequired": true, "weight": 0.1 }, {"filter": "Magnetizing Inductance", "invert": true, "log": false, "strictlyRequired": true, "weight": 1 }, {"filter": "Dc Current Density", "invert": true, "log": false, "strictlyRequired": false, "weight": 1 }, {"filter": "Effective Current Density", "invert": true, "log": false, "strictlyRequired": false, "weight": 1 }, {"filter": "Losses No Proximity", "invert": true, "log": false, "strictlyRequired": false, "weight": 1 }])");
        std::vector<OpenMagnetics::MagneticFilterOperation> filterFlow = filterFlowJson.get<std::vector<OpenMagnetics::MagneticFilterOperation>>();

        MagneticAdviser magneticAdviser(false);
        auto masMagnetics = magneticAdviser.get_advised_magnetic(mas.get_inputs(), magneticsCache.get(), filterFlow, 1);
        REQUIRE(masMagnetics.size() > 0);
    }

    TEST_CASE("Test_CatalogueAdviser_Web_5", "[adviser][magnetic-adviser][catalogue][bug]") {

        SKIP("Test needs investigation");
        auto catalogue_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "magnetics.ndjson").string();

        std::ifstream ndjsonFile(catalogue_path);
        std::string jsonLine;
        json coresJson = json::array();
        while (std::getline(ndjsonFile, jsonLine)) {
            json jf = json::parse(jsonLine);
            OpenMagnetics::Magnetic magnetic(jf); 
            magneticsCache.load(magnetic.get_reference(), magnetic);
        }
        magneticsCache.autocomplete_magnetics();

        auto json_path_2310 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_catalogueadviser_web_5_2310.json");
        std::ifstream json_file_2310(json_path_2310);
        json inputsJson = json::parse(json_file_2310);
        json filterFlowJson = json::parse(R"([{"filter":"Solid Insulation Requirements","invert":true,"log":false,"strictlyRequired":true,"weight":1},{"filter":"Turns Ratios","invert":true,"log":false,"strictlyRequired":true,"weight":0.1},{"filter":"Magnetizing Inductance","invert":true,"log":false,"strictlyRequired":false,"weight":1},{"filter":"Dc Current Density","invert":true,"log":false,"strictlyRequired":false,"weight":0.1},{"filter":"Effective Current Density","invert":true,"log":false,"strictlyRequired":false,"weight":0.1},{"filter":"Volume","invert":true,"log":false,"strictlyRequired":false,"weight":0.1},{"filter":"Area","invert":true,"log":false,"strictlyRequired":false,"weight":0.1},{"filter":"Height","invert":true,"log":false,"strictlyRequired":false,"weight":0.1},{"filter":"Losses No Proximity","invert":true,"log":false,"strictlyRequired":false,"weight":0.1}])");
        std::vector<OpenMagnetics::MagneticFilterOperation> filterFlow = filterFlowJson.get<std::vector<OpenMagnetics::MagneticFilterOperation>>();
        OpenMagnetics::Inputs inputs(inputsJson);

        MagneticAdviser magneticAdviser(false);
        auto masMagnetics = magneticAdviser.get_advised_magnetic(inputs, magneticsCache.get(), filterFlow, 9);
    }

    TEST_CASE("Test_CatalogueAdviser_Web_6", "[adviser][magnetic-adviser][catalogue][bug]") {

        auto catalogue_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "magnetics.ndjson").string();

        std::ifstream ndjsonFile(catalogue_path);
        std::string jsonLine;
        json coresJson = json::array();
        while (std::getline(ndjsonFile, jsonLine)) {
            json jf = json::parse(jsonLine);
            OpenMagnetics::Magnetic magnetic(jf); 
            magneticsCache.load(magnetic.get_reference(), magnetic);
        }
        magneticsCache.autocomplete_magnetics();

        json inputsJson = json::parse(R"({"designRequirements":{"application":null,"insulation":null,"isolationSides":["primary"],"leakageInductance":null,"magnetizingInductance":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":3.64583e-05,"nominal":null},"market":null,"maximumDimensions":null,"maximumWeight":null,"minimumImpedance":null,"name":null,"operatingTemperature":null,"strayCapacitance":null,"subApplication":null,"terminalType":null,"topology":"Buck Converter","turnsRatios":[],"wiringTechnology":null},"operatingPoints":[{"conditions":{"ambientRelativeHumidity":null,"ambientTemperature":25,"cooling":null,"name":null},"excitationsPerWinding":[{"current":{"harmonics":{"amplitudes":[1.999980766267635,0.34317551526427287,0.06653737321395338,0.015221897510571533,0.023296696518573023,0.007097693342598282,0.007156887875428565],"frequencies":[0,100000,200000,300000,400000,600000,700000]},"processed":{"acEffectiveFrequency":114080.4656948053,"average":2,"deadTime":null,"dutyCycle":0.6289062499999998,"effectiveFrequency":19716.313969852974,"label":"Triangular","negativePeak":1.570250137204429,"offset":1.9988219576129511,"peak":2.4273937780214734,"peakToPeak":0.8571436408170443,"phase":null,"positivePeak":2.4273937780214734,"rms":2.0153164343665915,"thd":0.21383578507534315},"waveform":{"ancillaryLabel":"Triangular","data":[1.570250137204429,2.4297498627955707,1.570250137204429],"numberPeriods":null,"time":[0,6.267179769103905e-06,1e-05]}},"frequency":100000,"magneticFieldStrength":null,"magneticFluxDensity":null,"magnetizingCurrent":null,"name":"Primary","voltage":{"harmonics":{"amplitudes":[0.018382352941177627,6.293943646607674,2.4093135980069746,0.869710676144456,1.705696569374191,0.8056909474220418,0.9034827028513644,0.7049759701548999,0.48653933430575286,0.575943248473444,0.3509141453482448,0.4291837721043859,0.38115834690001305,0.35466410125302295,0.28871363834165864,0.26825097600885595,0.2635392022165127,0.21628103381099675,0.17471047437045822,0.1679964488191641],"frequencies":[0,100000,200000,300000,400000,600000,700000,900000,1000000,1200000,1400000,1500000,1700000,2000000,2300000,2500000,2800000,3600000,5200000,6000000]},"processed":{"acEffectiveFrequency":619259.4938397474,"average":3.388131789017201e-16,"deadTime":null,"dutyCycle":0.6249999999999998,"effectiveFrequency":619257.5442914828,"label":"Rectangular","negativePeak":-6.7058823529411775,"offset":0,"peak":6.7058823529411775,"peakToPeak":10.7,"phase":null,"positivePeak":3.994117647058822,"rms":5.180147841606425,"thd":0.5956193445417157},"waveform":{"ancillaryLabel":"Rectangular","data":[-6.7058823529411775,3.994117647058822,3.994117647058822,-6.7058823529411775,-6.7058823529411775],"numberPeriods":null,"time":[0,0,6.267179769103905e-06,6.267179769103905e-06,1e-05]}}}],"name":"Min. input volt."},{"conditions":{"ambientRelativeHumidity":null,"ambientTemperature":25,"cooling":null,"name":null},"excitationsPerWinding":[{"current":{"harmonics":{"amplitudes":[1.9999699000112172,0.41060006607570226,0.018053537646037173,0.04423082015289921,0.008893668342075958,0.014935480811742638,0.0057827883819942735,0.006882280118733474,0.004185899621487101,0.0035882552214387634,0.00319653294622734,0.0019421043605148025,0.0025129518363500294,0.0010174460255423041,0.0020061769488597474,0.000460676521186724,0.001611989161582202,0.00012165486586560911,0.0012949290794106267,0.00013739575059009267,0.0010338975278241394,0.0002737513982975764,0.0008156147040044481,0.00036078251612133103,0.0006313615344976928,0.00040994461766320767,0.00047525425601085714,0.00043140362281954536,0.00034333614918288975,0.00043246933884406413,0.00023327234354635854,0.00041847872136654475,0.0001450732547597747,0.0003934595590190821,8.561203535289664e-05,0.0003605642211832911,7.533495403557772e-05,0.0003223566584951236,0.00010298615811853537,0.0002810151099167663,0.00013510958264978702,0.00023849746927257325,0.0001612655902563043,0.00019671622165927762,0.00017953354384482993,0.00015778907154269907,0.000189929520696439,0.00012445282619479882,0.0001930426019001372,0.0001005306008181055,0.00018970342158044007,9.009033812603397e-05,0.00018090272192525556,9.33309004608177e-05,0.00016779942203171963,0.0001048911519309436,0.00015178879849072472,0.00011877870872330486,0.0001346391276824486,0.00013125165787681355,0.0001186929371331471,0.0001403516089592326,0.00010697064153513654,0.00014510517561448124],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]},"processed":{"acEffectiveFrequency":110355.84175819751,"average":2,"deadTime":null,"dutyCycle":0.5273437499999998,"effectiveFrequency":22365.154374948852,"label":"Triangular","negativePeak":1.4930981933324685,"offset":1.9955986527614609,"peak":2.4999685937819276,"peakToPeak":1.006870400449459,"phase":null,"positivePeak":2.4999685937819276,"rms":2.0212720866579925,"thd":0.127253075838145},"waveform":{"ancillaryLabel":"Triangular","data":[1.4930981933324685,2.5069018066675315,1.4930981933324685],"numberPeriods":null,"time":[0,5.280222325150533e-06,1e-05]}},"frequency":100000,"magneticFieldStrength":null,"magneticFluxDensity":null,"magnetizingCurrent":null,"name":"Primary","voltage":{"harmonics":{"amplitudes":[0.058226102941177,8.06397384761784,0.5934014383980282,2.631917895729841,0.5876866581600322,1.512455986518031,0.5782231957437205,1.0099203310757312,0.5651021894654171,0.713996388702092,0.5484500017882536,0.5130264601577592,0.5284270023825729,0.36426072438637475,0.5052260236783366,0.24790491440874138,0.4790705037834115,0.15360284899194315,0.4502123346524472,0.07545461504925513,0.41892943622965506,0.009880133690457456,0.38552307992781826,0.04540369425376326,0.3503149872198955,0.09191569970740314,0.3136442312853652,0.130706996323082,0.27586397155015135,0.16253521765579576,0.23733805256829477,0.18797152250570992,0.1984375,0.2074687410079629,0.1595369474317056,0.22140606271007607,0.12101102844984904,0.23011900864900559,0.0832307687146352,0.2339198453541049,0.0465600127801049,0.23311158968563056,0.011351920072182047,0.22799758686455546,0.022054436229654827,0.2188879455024585,0.05333733465244702,0.2061036837125316,0.08219550378341134,0.18997916971048465,0.10835102367833661,0.17086326627866774,0.1315520023825731,0.1491194745581959,0.15157500178825387,0.12512529695002444,0.16822718946541748,0.09927098800643068,0.181348195743721,0.07195782766572233,0.19081165816003276,0.04359602766923337,0.19652643839802877,0.014602366048946653],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]},"processed":{"acEffectiveFrequency":592064.6087479359,"average":-6.776263578034402e-16,"deadTime":null,"dutyCycle":0.5234374999999998,"effectiveFrequency":592052.1336211944,"label":"Rectangular","negativePeak":-6.705882352941177,"offset":0,"peak":6.705882352941177,"peakToPeak":12.7,"phase":null,"positivePeak":5.994117647058823,"rms":6.343287084461784,"thd":0.48696668674418153},"waveform":{"ancillaryLabel":"Rectangular","data":[-6.705882352941177,5.994117647058823,5.994117647058823,-6.705882352941177,-6.705882352941177],"numberPeriods":null,"time":[0,0,5.280222325150533e-06,5.280222325150533e-06,1e-05]}}}],"name":"Max. input volt."}]})");
        json filterFlowJson = json::parse(R"([{"filter":"Turns Ratios","invert":true,"log":false,"strictlyRequired":true,"weight":1},{"filter":"Solid Insulation Requirements","invert":true,"log":false,"strictlyRequired":true,"weight":1},{"filter":"Magnetizing Inductance","invert":true,"log":false,"strictlyRequired":true,"weight":1},{"filter":"Dc Current Density","invert":true,"log":false,"strictlyRequired":false,"weight":0.1},{"filter":"Effective Current Density","invert":true,"log":false,"strictlyRequired":false,"weight":0.1},{"filter":"Volume","invert":true,"log":false,"strictlyRequired":false,"weight":0.1},{"filter":"Area","invert":true,"log":false,"strictlyRequired":false,"weight":0.1},{"filter":"Height","invert":true,"log":false,"strictlyRequired":false,"weight":0.1},{"filter":"Losses No Proximity","invert":true,"log":false,"strictlyRequired":false,"weight":0.1}])");
        std::vector<OpenMagnetics::MagneticFilterOperation> filterFlow = filterFlowJson.get<std::vector<OpenMagnetics::MagneticFilterOperation>>();
        OpenMagnetics::Inputs inputs(inputsJson);

        MagneticAdviser magneticAdviser(false);
        auto masMagnetics = magneticAdviser.get_advised_magnetic(inputs, magneticsCache.get(), filterFlow, 9);
    }

    TEST_CASE("Test_CatalogueAdviser_Web_7", "[adviser][magnetic-adviser][catalogue][bug]") {

        auto catalogue_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "magnetics.ndjson").string();

        std::ifstream ndjsonFile(catalogue_path);
        std::string jsonLine;
        json coresJson = json::array();
        while (std::getline(ndjsonFile, jsonLine)) {
            json jf = json::parse(jsonLine);
            OpenMagnetics::Magnetic magnetic(jf); 
            magneticsCache.load(magnetic.get_reference(), magnetic);
        }
        magneticsCache.autocomplete_magnetics();

        json inputsJson = json::parse(R"({"designRequirements":{"application":null,"insulation":null,"isolationSides":["primary"],"leakageInductance":null,"magnetizingInductance":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":3.64583e-05,"nominal":null},"market":null,"maximumDimensions":null,"maximumWeight":null,"minimumImpedance":null,"name":null,"operatingTemperature":null,"strayCapacitance":null,"subApplication":null,"terminalType":null,"topology":"Buck Converter","turnsRatios":[],"wiringTechnology":null},"operatingPoints":[{"conditions":{"ambientRelativeHumidity":null,"ambientTemperature":25,"cooling":null,"name":null},"excitationsPerWinding":[{"current":{"harmonics":{"amplitudes":[1.999980766267635,0.34317551526427287,0.06653737321395338,0.015221897510571533,0.023296696518573023,0.007097693342598282,0.007156887875428565],"frequencies":[0,100000,200000,300000,400000,600000,700000]},"processed":{"acEffectiveFrequency":114080.4656948053,"average":2,"deadTime":null,"dutyCycle":0.6289062499999998,"effectiveFrequency":19716.313969852974,"label":"Triangular","negativePeak":1.570250137204429,"offset":1.9988219576129511,"peak":2.4273937780214734,"peakToPeak":0.8571436408170443,"phase":null,"positivePeak":2.4273937780214734,"rms":2.0153164343665915,"thd":0.21383578507534315},"waveform":{"ancillaryLabel":"Triangular","data":[1.570250137204429,2.4297498627955707,1.570250137204429],"numberPeriods":null,"time":[0,6.267179769103905e-06,1e-05]}},"frequency":100000,"magneticFieldStrength":null,"magneticFluxDensity":null,"magnetizingCurrent":null,"name":"Primary","voltage":{"harmonics":{"amplitudes":[0.018382352941177627,6.293943646607674,2.4093135980069746,0.869710676144456,1.705696569374191,0.8056909474220418,0.9034827028513644,0.7049759701548999,0.48653933430575286,0.575943248473444,0.3509141453482448,0.4291837721043859,0.38115834690001305,0.35466410125302295,0.28871363834165864,0.26825097600885595,0.2635392022165127,0.21628103381099675,0.17471047437045822,0.1679964488191641],"frequencies":[0,100000,200000,300000,400000,600000,700000,900000,1000000,1200000,1400000,1500000,1700000,2000000,2300000,2500000,2800000,3600000,5200000,6000000]},"processed":{"acEffectiveFrequency":619259.4938397474,"average":3.388131789017201e-16,"deadTime":null,"dutyCycle":0.6249999999999998,"effectiveFrequency":619257.5442914828,"label":"Rectangular","negativePeak":-6.7058823529411775,"offset":0,"peak":6.7058823529411775,"peakToPeak":10.7,"phase":null,"positivePeak":3.994117647058822,"rms":5.180147841606425,"thd":0.5956193445417157},"waveform":{"ancillaryLabel":"Rectangular","data":[-6.7058823529411775,3.994117647058822,3.994117647058822,-6.7058823529411775,-6.7058823529411775],"numberPeriods":null,"time":[0,0,6.267179769103905e-06,6.267179769103905e-06,1e-05]}}}],"name":"Min. input volt."},{"conditions":{"ambientRelativeHumidity":null,"ambientTemperature":25,"cooling":null,"name":null},"excitationsPerWinding":[{"current":{"harmonics":{"amplitudes":[1.9999699000112172,0.41060006607570226,0.018053537646037173,0.04423082015289921,0.008893668342075958,0.014935480811742638,0.0057827883819942735,0.006882280118733474,0.004185899621487101,0.0035882552214387634,0.00319653294622734,0.0019421043605148025,0.0025129518363500294,0.0010174460255423041,0.0020061769488597474,0.000460676521186724,0.001611989161582202,0.00012165486586560911,0.0012949290794106267,0.00013739575059009267,0.0010338975278241394,0.0002737513982975764,0.0008156147040044481,0.00036078251612133103,0.0006313615344976928,0.00040994461766320767,0.00047525425601085714,0.00043140362281954536,0.00034333614918288975,0.00043246933884406413,0.00023327234354635854,0.00041847872136654475,0.0001450732547597747,0.0003934595590190821,8.561203535289664e-05,0.0003605642211832911,7.533495403557772e-05,0.0003223566584951236,0.00010298615811853537,0.0002810151099167663,0.00013510958264978702,0.00023849746927257325,0.0001612655902563043,0.00019671622165927762,0.00017953354384482993,0.00015778907154269907,0.000189929520696439,0.00012445282619479882,0.0001930426019001372,0.0001005306008181055,0.00018970342158044007,9.009033812603397e-05,0.00018090272192525556,9.33309004608177e-05,0.00016779942203171963,0.0001048911519309436,0.00015178879849072472,0.00011877870872330486,0.0001346391276824486,0.00013125165787681355,0.0001186929371331471,0.0001403516089592326,0.00010697064153513654,0.00014510517561448124],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]},"processed":{"acEffectiveFrequency":110355.84175819751,"average":2,"deadTime":null,"dutyCycle":0.5273437499999998,"effectiveFrequency":22365.154374948852,"label":"Triangular","negativePeak":1.4930981933324685,"offset":1.9955986527614609,"peak":2.4999685937819276,"peakToPeak":1.006870400449459,"phase":null,"positivePeak":2.4999685937819276,"rms":2.0212720866579925,"thd":0.127253075838145},"waveform":{"ancillaryLabel":"Triangular","data":[1.4930981933324685,2.5069018066675315,1.4930981933324685],"numberPeriods":null,"time":[0,5.280222325150533e-06,1e-05]}},"frequency":100000,"magneticFieldStrength":null,"magneticFluxDensity":null,"magnetizingCurrent":null,"name":"Primary","voltage":{"harmonics":{"amplitudes":[0.058226102941177,8.06397384761784,0.5934014383980282,2.631917895729841,0.5876866581600322,1.512455986518031,0.5782231957437205,1.0099203310757312,0.5651021894654171,0.713996388702092,0.5484500017882536,0.5130264601577592,0.5284270023825729,0.36426072438637475,0.5052260236783366,0.24790491440874138,0.4790705037834115,0.15360284899194315,0.4502123346524472,0.07545461504925513,0.41892943622965506,0.009880133690457456,0.38552307992781826,0.04540369425376326,0.3503149872198955,0.09191569970740314,0.3136442312853652,0.130706996323082,0.27586397155015135,0.16253521765579576,0.23733805256829477,0.18797152250570992,0.1984375,0.2074687410079629,0.1595369474317056,0.22140606271007607,0.12101102844984904,0.23011900864900559,0.0832307687146352,0.2339198453541049,0.0465600127801049,0.23311158968563056,0.011351920072182047,0.22799758686455546,0.022054436229654827,0.2188879455024585,0.05333733465244702,0.2061036837125316,0.08219550378341134,0.18997916971048465,0.10835102367833661,0.17086326627866774,0.1315520023825731,0.1491194745581959,0.15157500178825387,0.12512529695002444,0.16822718946541748,0.09927098800643068,0.181348195743721,0.07195782766572233,0.19081165816003276,0.04359602766923337,0.19652643839802877,0.014602366048946653],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]},"processed":{"acEffectiveFrequency":592064.6087479359,"average":-6.776263578034402e-16,"deadTime":null,"dutyCycle":0.5234374999999998,"effectiveFrequency":592052.1336211944,"label":"Rectangular","negativePeak":-6.705882352941177,"offset":0,"peak":6.705882352941177,"peakToPeak":12.7,"phase":null,"positivePeak":5.994117647058823,"rms":6.343287084461784,"thd":0.48696668674418153},"waveform":{"ancillaryLabel":"Rectangular","data":[-6.705882352941177,5.994117647058823,5.994117647058823,-6.705882352941177,-6.705882352941177],"numberPeriods":null,"time":[0,0,5.280222325150533e-06,5.280222325150533e-06,1e-05]}}}],"name":"Max. input volt."}]})");
        json filterFlowJson = json::parse(R"([{"filter":"Turns Ratios","invert":true,"log":false,"strictlyRequired":true,"weight":1},{"filter":"Solid Insulation Requirements","invert":true,"log":false,"strictlyRequired":true,"weight":1},{"filter":"Magnetizing Inductance","invert":true,"log":false,"strictlyRequired":true,"weight":1},{"filter":"Dc Current Density","invert":true,"log":false,"strictlyRequired":false,"weight":0.1},{"filter":"Effective Current Density","invert":true,"log":false,"strictlyRequired":false,"weight":0.1},{"filter":"Volume","invert":true,"log":false,"strictlyRequired":false,"weight":0.1},{"filter":"Area","invert":true,"log":false,"strictlyRequired":false,"weight":0.1},{"filter":"Height","invert":true,"log":false,"strictlyRequired":false,"weight":0.1},{"filter":"Losses No Proximity","invert":true,"log":false,"strictlyRequired":false,"weight":0.1}])");
        std::vector<OpenMagnetics::MagneticFilterOperation> filterFlow = filterFlowJson.get<std::vector<OpenMagnetics::MagneticFilterOperation>>();
        OpenMagnetics::Inputs inputs(inputsJson);

        MagneticAdviser magneticAdviser(false);
        auto masMagnetics = magneticAdviser.get_advised_magnetic(inputs, magneticsCache.get(), filterFlow, 9);
    }

    TEST_CASE("Test_CatalogueAdviser_Web_8", "[adviser][magnetic-adviser][catalogue][bug]") {

        auto catalogue_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "magnetics.ndjson").string();

        std::ifstream ndjsonFile(catalogue_path);
        std::string jsonLine;
        json coresJson = json::array();
        while (std::getline(ndjsonFile, jsonLine)) {
            json jf = json::parse(jsonLine);
            OpenMagnetics::Magnetic magnetic(jf); 
            magneticsCache.load(magnetic.get_reference(), magnetic);
        }
        magneticsCache.autocomplete_magnetics();

        json inputsJson = json::parse(R"({"designRequirements":{"application":null,"insulation":null,"isolationSides":["primary"],"leakageInductance":null,"magnetizingInductance":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":3.64583e-05,"nominal":null},"market":null,"maximumDimensions":null,"maximumWeight":null,"minimumImpedance":null,"name":null,"operatingTemperature":null,"strayCapacitance":null,"subApplication":null,"terminalType":null,"topology":"Buck Converter","turnsRatios":[],"wiringTechnology":null},"operatingPoints":[{"conditions":{"ambientRelativeHumidity":null,"ambientTemperature":25,"cooling":null,"name":null},"excitationsPerWinding":[{"current":{"harmonics":{"amplitudes":[1.999980766267635,0.34317551526427287,0.06653737321395338,0.015221897510571533,0.023296696518573023,0.006080616319137073,0.007097693342598282,0.007156887875428565,0.00025275988102715387,0.004190061229987666,0.002795472797738478,0.0010230934707582443,0.0026168865956771228,0.0009964611944203304,0.0012663235380700245,0.001608689644721802,0.00012885586992262952,0.0011764174617886498,0.0009127330463713986,0.00031522965171201,0.0009642338598537787,0.0004221986598760322,0.0005074249547744154,0.0007138584136362339,8.875746692601736e-05,0.0005543955849936191,0.0004662280350403915,0.00015098833800040892,0.0005114159601172801,0.0002445643993118587,0.0002782714950257878,0.00041538128182986206,6.973629428293869e-05,0.0003323457336371823,0.0002931857972092662,9.606860343299165e-05,0.0003284632220442052,0.00016608689599784698,0.00018527071146608962,0.00028230144044252636,5.9305843355891194e-05,0.00023295463002489185,0.00020913069335624248,7.701851433203897e-05,0.00024044614725589024,0.00012445718877296657,0.000143851803955198,0.00021427758105438605,5.33738489132515e-05,0.0001851670688564565,0.0001634765278633521,7.279356736110384e-05,0.00019613657883187313,9.997823242314758e-05,0.00012761422865876747,0.00017868499827029954,5.027706631361615e-05,0.00016494308009893058,0.00013827984056522845,7.61005492886236e-05,0.0001771726454490597,8.492772813281776e-05,0.0001267766955395071,0.00016343008633114346],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]},"processed":{"acEffectiveFrequency":114080.4656948053,"average":2,"deadTime":null,"dutyCycle":0.6289062499999998,"effectiveFrequency":19716.313969852974,"label":"Triangular","negativePeak":1.570250137204429,"offset":1.9988219576129511,"peak":2.4273937780214734,"peakToPeak":0.8571436408170443,"phase":null,"positivePeak":2.4273937780214734,"rms":2.0153164343665915,"thd":0.21383578507534315},"waveform":{"ancillaryLabel":"Triangular","data":[1.570250137204429,2.4297498627955707,1.570250137204429],"numberPeriods":null,"time":[0,6.267179769103905e-06,1e-05]}},"frequency":100000,"magneticFieldStrength":null,"magneticFluxDensity":null,"magnetizingCurrent":null,"name":"Primary","voltage":{"harmonics":{"amplitudes":[0.018382352941177627,6.293943646607674,2.4093135980069746,0.869710676144456,1.705696569374191,0.5226659051133817,0.8056909474220418,0.9034827028513644,0,0.7049759701548999,0.48653933430575286,0.23988311228254675,0.575943248473444,0.20396433106746678,0.3509141453482448,0.4291837721043859,0,0.38115834690001305,0.27650100956849705,0.14230043179225071,0.35466410125302295,0.1298034510110597,0.22995289621561393,0.28871363834165864,0,0.26825097600885595,0.1984548480887606,0.10399317481106156,0.2635392022165127,0.09795245935089585,0.17603728477800934,0.2240058405856783,0,0.21327128956083283,0.15955089470036488,0.08449437251521935,0.21628103381099675,0.08115707017428392,0.14718414649536288,0.18892365276003928,0,0.18282589515507067,0.13782846182499203,0.07353274672801598,0.1895719116417196,0.07162801802739381,0.13077523802428803,0.1689558751458146,0,0.16555444965346452,0.125559088895798,0.06738071355457714,0.17471047437045822,0.06638459764585816,0.1218717588637492,0.15830764798641112,0,0.15676912039455004,0.119512962336555,0.0644646908254962,0.1679964488191641,0.06415371361480038,0.11836198717682256,0.15450764415544221],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]},"processed":{"acEffectiveFrequency":619259.4938397474,"average":3.388131789017201e-16,"deadTime":null,"dutyCycle":0.6249999999999998,"effectiveFrequency":619257.5442914828,"label":"Rectangular","negativePeak":-6.7058823529411775,"offset":0,"peak":6.7058823529411775,"peakToPeak":10.7,"phase":null,"positivePeak":3.994117647058822,"rms":5.180147841606425,"thd":0.5956193445417157},"waveform":{"ancillaryLabel":"Rectangular","data":[-6.7058823529411775,3.994117647058822,3.994117647058822,-6.7058823529411775,-6.7058823529411775],"numberPeriods":null,"time":[0,0,6.267179769103905e-06,6.267179769103905e-06,1e-05]}}}],"name":"Min. input volt."},{"conditions":{"ambientRelativeHumidity":null,"ambientTemperature":25,"cooling":null,"name":null},"excitationsPerWinding":[{"current":{"harmonics":{"amplitudes":[1.9999699000112172,0.41060006607570226,0.018053537646037173,0.04423082015289921,0.008893668342075958,0.014935480811742638,0.0057827883819942735,0.006882280118733474,0.004185899621487101,0.0035882552214387634,0.00319653294622734,0.0019421043605148025,0.0025129518363500294,0.0010174460255423041,0.0020061769488597474,0.000460676521186724,0.001611989161582202,0.00012165486586560911,0.0012949290794106267,0.00013739575059009267,0.0010338975278241394,0.0002737513982975764,0.0008156147040044481,0.00036078251612133103,0.0006313615344976928,0.00040994461766320767,0.00047525425601085714,0.00043140362281954536,0.00034333614918288975,0.00043246933884406413,0.00023327234354635854,0.00041847872136654475,0.0001450732547597747,0.0003934595590190821,8.561203535289664e-05,0.0003605642211832911,7.533495403557772e-05,0.0003223566584951236,0.00010298615811853537,0.0002810151099167663,0.00013510958264978702,0.00023849746927257325,0.0001612655902563043,0.00019671622165927762,0.00017953354384482993,0.00015778907154269907,0.000189929520696439,0.00012445282619479882,0.0001930426019001372,0.0001005306008181055,0.00018970342158044007,9.009033812603397e-05,0.00018090272192525556,9.33309004608177e-05,0.00016779942203171963,0.0001048911519309436,0.00015178879849072472,0.00011877870872330486,0.0001346391276824486,0.00013125165787681355,0.0001186929371331471,0.0001403516089592326,0.00010697064153513654,0.00014510517561448124],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]},"processed":{"acEffectiveFrequency":110355.84175819751,"average":2,"deadTime":null,"dutyCycle":0.5273437499999998,"effectiveFrequency":22365.154374948852,"label":"Triangular","negativePeak":1.4930981933324685,"offset":1.9955986527614609,"peak":2.4999685937819276,"peakToPeak":1.006870400449459,"phase":null,"positivePeak":2.4999685937819276,"rms":2.0212720866579925,"thd":0.127253075838145},"waveform":{"ancillaryLabel":"Triangular","data":[1.4930981933324685,2.5069018066675315,1.4930981933324685],"numberPeriods":null,"time":[0,5.280222325150533e-06,1e-05]}},"frequency":100000,"magneticFieldStrength":null,"magneticFluxDensity":null,"magnetizingCurrent":null,"name":"Primary","voltage":{"harmonics":{"amplitudes":[0.058226102941177,8.06397384761784,0.5934014383980282,2.631917895729841,0.5876866581600322,1.512455986518031,0.5782231957437205,1.0099203310757312,0.5651021894654171,0.713996388702092,0.5484500017882536,0.5130264601577592,0.5284270023825729,0.36426072438637475,0.5052260236783366,0.24790491440874138,0.4790705037834115,0.15360284899194315,0.4502123346524472,0.07545461504925513,0.41892943622965506,0.009880133690457456,0.38552307992781826,0.04540369425376326,0.3503149872198955,0.09191569970740314,0.3136442312853652,0.130706996323082,0.27586397155015135,0.16253521765579576,0.23733805256829477,0.18797152250570992,0.1984375,0.2074687410079629,0.1595369474317056,0.22140606271007607,0.12101102844984904,0.23011900864900559,0.0832307687146352,0.2339198453541049,0.0465600127801049,0.23311158968563056,0.011351920072182047,0.22799758686455546,0.022054436229654827,0.2188879455024585,0.05333733465244702,0.2061036837125316,0.08219550378341134,0.18997916971048465,0.10835102367833661,0.17086326627866774,0.1315520023825731,0.1491194745581959,0.15157500178825387,0.12512529695002444,0.16822718946541748,0.09927098800643068,0.181348195743721,0.07195782766572233,0.19081165816003276,0.04359602766923337,0.19652643839802877,0.014602366048946653],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]},"processed":{"acEffectiveFrequency":592064.6087479359,"average":-6.776263578034402e-16,"deadTime":null,"dutyCycle":0.5234374999999998,"effectiveFrequency":592052.1336211944,"label":"Rectangular","negativePeak":-6.705882352941177,"offset":0,"peak":6.705882352941177,"peakToPeak":12.7,"phase":null,"positivePeak":5.994117647058823,"rms":6.343287084461784,"thd":0.48696668674418153},"waveform":{"ancillaryLabel":"Rectangular","data":[-6.705882352941177,5.994117647058823,5.994117647058823,-6.705882352941177,-6.705882352941177],"numberPeriods":null,"time":[0,0,5.280222325150533e-06,5.280222325150533e-06,1e-05]}}}],"name":"Max. input volt."}]})");
        json filterFlowJson = json::parse(R"([{"filter":"Turns Ratios","invert":true,"log":false,"strictlyRequired":true,"weight":1},{"filter":"Solid Insulation Requirements","invert":true,"log":false,"strictlyRequired":true,"weight":1},{"filter":"Magnetizing Inductance","invert":true,"log":false,"strictlyRequired":true,"weight":1},{"filter":"Dc Current Density","invert":true,"log":false,"strictlyRequired":false,"weight":0.1},{"filter":"Effective Current Density","invert":true,"log":false,"strictlyRequired":false,"weight":0.1},{"filter":"Volume","invert":true,"log":false,"strictlyRequired":false,"weight":0.1},{"filter":"Area","invert":true,"log":false,"strictlyRequired":false,"weight":0.1},{"filter":"Height","invert":true,"log":false,"strictlyRequired":false,"weight":0.1},{"filter":"Losses No Proximity","invert":true,"log":false,"strictlyRequired":false,"weight":0.1}])");
        std::vector<OpenMagnetics::MagneticFilterOperation> filterFlow = filterFlowJson.get<std::vector<OpenMagnetics::MagneticFilterOperation>>();
        OpenMagnetics::Inputs inputs(inputsJson);

        MagneticAdviser magneticAdviser(false);
        auto masMagnetics = magneticAdviser.get_advised_magnetic(inputs, magneticsCache.get(), filterFlow, 9);
    }

    TEST_CASE("Test_CatalogueAdviser_Power_Found", "[adviser][magnetic-adviser][catalogue]") { 
        json isolatedbuckInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 36;
        inputVoltage["maximum"] = 72;
        isolatedbuckInputsJson["inputVoltage"] = inputVoltage;
        isolatedbuckInputsJson["diodeVoltageDrop"] = 0.7;
        isolatedbuckInputsJson["maximumSwitchCurrent"] = 0.7;
        isolatedbuckInputsJson["efficiency"] = 1;
        isolatedbuckInputsJson["operatingPoints"] = json::array();

        {
            json isolatedbuckOperatingPointJson;
            isolatedbuckOperatingPointJson["outputVoltages"] = {10, 10};
            isolatedbuckOperatingPointJson["outputCurrents"] = {0.02, 0.1};
            isolatedbuckOperatingPointJson["switchingFrequency"] = 750000;
            isolatedbuckOperatingPointJson["ambientTemperature"] = 42;
            isolatedbuckInputsJson["operatingPoints"].push_back(isolatedbuckOperatingPointJson);
        }
        OpenMagnetics::IsolatedBuck isolatedbuckInputs(isolatedbuckInputsJson);
        isolatedbuckInputs._assertErrors = true;

        auto inputs = isolatedbuckInputs.process();

        auto catalogue_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "magnetics.ndjson").string();

        std::ifstream ndjsonFile(catalogue_path);
        std::string jsonLine;
        json coresJson = json::array();
        while (std::getline(ndjsonFile, jsonLine)) {
            json jf = json::parse(jsonLine);
            OpenMagnetics::Magnetic magnetic(jf); 
            magneticsCache.load(magnetic.get_reference(), magnetic);
        }
        magneticsCache.autocomplete_magnetics();

        std::vector<OpenMagnetics::MagneticFilterOperation> flow{
            MagneticFilterOperation(OpenMagnetics::MagneticFilters::TURNS_RATIOS, true, false, true, 1.0),
            MagneticFilterOperation(OpenMagnetics::MagneticFilters::MAXIMUM_DIMENSIONS, true, false, 1.0),
            MagneticFilterOperation(OpenMagnetics::MagneticFilters::SATURATION, true, false, 1.0),
            MagneticFilterOperation(OpenMagnetics::MagneticFilters::DC_CURRENT_DENSITY, true, false, 1.0),
            MagneticFilterOperation(OpenMagnetics::MagneticFilters::EFFECTIVE_CURRENT_DENSITY, true, false, 1.0),
            MagneticFilterOperation(OpenMagnetics::MagneticFilters::MAGNETIZING_INDUCTANCE, true, false, 1.0),
        };

        MagneticAdviser magneticAdviser;
        auto masMagnetics = magneticAdviser.get_advised_magnetic(inputs, magneticsCache.get(), flow, 1);
        REQUIRE(masMagnetics.size() > 0);

        auto scorings = magneticAdviser.get_scorings();
        for (auto [name, values] : scorings) {
            double scoringTotal = 0;
            for (auto [key, value] : values) {
                scoringTotal += value;
            }
        }

        for (auto [mas, scoring] : masMagnetics) {

            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            {
                auto outFile = outputFilePath;
                std::string filename = "MagneticAdviser" +  std::to_string(scoring) + ".mas.json";
                outFile.append(filename);
                to_file(outFile, mas);
            }

            auto operatingPoint = inputs.get_operating_points()[0];
            MagnetizingInductance magnetizingInductanceObj;
            auto magneticFluxDensity = magnetizingInductanceObj.calculate_inductance_and_magnetic_flux_density(mas.get_mutable_magnetic().get_core(), mas.get_mutable_magnetic().get_coil(), &operatingPoint).second;

            OpenMagneticsTesting::check_turns_description(mas.get_mutable_magnetic().get_coil());
            if (plot) {
                auto outFile = outputFilePath;
                std::string filename = "MagneticAdviser" + mas.get_magnetic().get_manufacturer_info().value().get_reference().value() + ".svg";
                outFile.append(filename);
                Painter painter(outFile, true);

                painter.paint_magnetic_field(mas.get_mutable_inputs().get_operating_point(0), mas.get_mutable_magnetic());
                painter.paint_core(mas.get_mutable_magnetic());
                painter.paint_bobbin(mas.get_mutable_magnetic());
                painter.paint_coil_turns(mas.get_mutable_magnetic());
                #ifdef ENABLE_MATPLOTPP
                painter.export_svg();
                #else
                    INFO("matplotplusplus disabled — skipping AdvancedPainter call");
                #endif

            }
        }
    }

    TEST_CASE("Test_MagneticAdviser_Flyback_Standard_Cores_Performance", "[adviser][magnetic-adviser][performance]") {
        // This test investigates performance issues with standard cores mode
        // for a flyback converter design
        clear_databases();
        settings.reset();

        auto testStart = std::chrono::high_resolution_clock::now();
        auto printElapsed = [&testStart](const std::string& msg) {
            auto now = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - testStart);
            std::cout << "[" << elapsed.count() << " ms] " << msg << std::endl;
            std::cout.flush();
        };

        printElapsed("Test started");

        // JSON inputs from Python example - flyback with 500uH magnetizing inductance
        std::string inputsString = R"({
            "designRequirements": {
                "magnetizingInductance": {
                    "nominal": 500e-6
                },
                "turnsRatios": [
                    {
                        "nominal": 0.1
                    }
                ],
                "topology": "Flyback Converter"
            },
            "operatingPoints": [
                {
                    "conditions": {
                        "ambientTemperature": 25
                    },
                    "excitationsPerWinding": [
                        {
                            "current": {
                                "waveform": {
                                    "data": [
                                        0,
                                        0.21208333333333335,
                                        0.45458333333333334,
                                        0.1425,
                                        0,
                                        0,
                                        0
                                    ],
                                    "time": [
                                        0,
                                        0,
                                        0.0000021375,
                                        0.0000021375,
                                        0.000004275,
                                        0.000005,
                                        0.000005
                                    ]
                                }
                            },
                            "frequency": 200000,
                            "name": "First primary",
                            "voltage": {
                                "waveform": {
                                    "data": [
                                        0,
                                        100,
                                        100,
                                        -101.4,
                                        -101.4,
                                        0,
                                        0
                                    ],
                                    "time": [
                                        0,
                                        0,
                                        0.0000021375,
                                        0.0000021375,
                                        0.000004275,
                                        0.000004275,
                                        0.000005
                                    ]
                                }
                            }
                        },
                        {
                            "current": {
                                "waveform": {
                                    "data": [
                                        0,
                                        4.25,
                                        5.75,
                                        0,
                                        0
                                    ],
                                    "time": [
                                        0,
                                        0,
                                        0.00000225,
                                        0.00000225,
                                        0.000005
                                    ]
                                }
                            },
                            "frequency": 200000,
                            "name": "Secondary 0",
                            "voltage": {
                                "waveform": {
                                    "data": [
                                        0,
                                        7.291333333333334,
                                        7.291333333333334,
                                        -6.135333333333334,
                                        -6.135333333333334,
                                        0,
                                        0
                                    ],
                                    "time": [
                                        0,
                                        0,
                                        0.00000225,
                                        0.00000225,
                                        0.000004275,
                                        0.000004275,
                                        0.000005
                                    ]
                                }
                            }
                        }
                    ]
                }
            ]
        })";

        printElapsed("Parsing JSON inputs");
        json inputsJson = json::parse(inputsString);
        OpenMagnetics::Inputs inputs(inputsJson);
        
        printElapsed("Processing inputs");
        inputs.process();
        printElapsed("Inputs processed");

        // Now test the full MagneticAdviser
        printElapsed("Creating MagneticAdviser");
        MagneticAdviser magneticAdviser;
        magneticAdviser.set_core_mode(CoreAdviser::CoreAdviserModes::STANDARD_CORES);
        
        printElapsed("Calling get_advised_magnetic with 5 results");
        auto masMagnetics = magneticAdviser.get_advised_magnetic(inputs, 5);
        printElapsed("MagneticAdviser returned " + std::to_string(masMagnetics.size()) + " magnetics");

        for (auto& [mas, scoring] : masMagnetics) {
            std::cout << "  Magnetic: " << mas.get_magnetic().get_core().get_name().value() << " - Score: " << scoring << std::endl;
        }

        printElapsed("Test completed");
    }

}  // namespace
