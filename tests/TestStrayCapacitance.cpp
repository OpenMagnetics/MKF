#include "StrayCapacitance.h"
#include "Settings.h"
#include "MasWrapper.h"
#include "Utils.h"
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


SUITE(StrayCapacitance) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    double maximumError = 0.2;
    auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
    bool plot = false;

    TEST(Test_Get_Surrounding_Turns) {
        settings->reset();
        std::string file_path = __FILE__;
        auto inventory_path = file_path.substr(0, file_path.rfind("/")).append("/testData/dont_worry_this_is_just_a_dummy_mas.json");
        std::ifstream jsonFile(inventory_path);
        json masJson;
        jsonFile >> masJson;

        auto mas = OpenMagnetics::MasWrapper(masJson);

        auto turns = mas.get_magnetic().get_coil().get_turns_description().value();

        auto surroundingTurns = OpenMagnetics::StrayCapacitance::get_surrounding_turns(turns[112], turns);
        std::cout << surroundingTurns.size() << std::endl;
        std::vector<OpenMagnetics::Turn> newTurns = surroundingTurns;
        newTurns.push_back(turns[112]);

        mas.get_mutable_magnetic().get_mutable_coil().set_turns_description(newTurns);

        auto outFile = outputFilePath;
        outFile.append("Test_Get_Surrounding_Turns.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        painter.paint_core(mas.get_mutable_magnetic());
        painter.paint_bobbin(mas.get_mutable_magnetic());
        painter.paint_coil_turns(mas.get_mutable_magnetic());
        painter.export_svg();
    }

    TEST(Test_Get_Voltage_Per_Turn) {
        settings->reset();
        std::string file_path = __FILE__;
        auto inventory_path = file_path.substr(0, file_path.rfind("/")).append("/testData/dont_worry_this_is_just_a_dummy_mas.json");
        std::ifstream jsonFile(inventory_path);
        json masJson;
        jsonFile >> masJson;

        auto mas = OpenMagnetics::MasWrapper(masJson);

        auto coil = mas.get_mutable_magnetic().get_mutable_coil();
        auto operatingPoint = mas.get_mutable_inputs().get_operating_point(0);

        auto strayCapacitanceOutput = OpenMagnetics::StrayCapacitance::calculate_voltages_per_turn(coil, operatingPoint);
        auto voltageDividerEndPerTurn = strayCapacitanceOutput.get_voltage_divider_end_per_turn().value();
        auto voltageDividerStartPerTurn = strayCapacitanceOutput.get_voltage_divider_start_per_turn().value();
        auto voltagePerTurn = strayCapacitanceOutput.get_voltage_per_turn().value();
        CHECK(voltageDividerEndPerTurn.size() == coil.get_turns_description().value().size());
        CHECK(voltageDividerStartPerTurn.size() == coil.get_turns_description().value().size());
        CHECK(voltagePerTurn.size() == coil.get_turns_description().value().size());
    }


    // TEST(Test_Get_Capacitance_Among_Turns) {
    //     settings->reset();
    //     json masJson = json::parse(R"({"outputs": [], "inputs": {"designRequirements": {"isolationSides": ["primary", "secondary"], "magnetizingInductance": {"nominal": 0.00039999999999999996 }, "name": "My Design Requirements", "turnsRatios": [] }, "operatingPoints": [{"conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"frequency": 100000, "current": {"processed": {"label": "Triangular", "peakToPeak": 0.5, "offset": 0, "dutyCycle": 0.5 } }, "voltage": {"processed": {"label": "Rectangular", "peakToPeak": 1000, "offset": 0, "dutyCycle": 0.5 } } } ], "name": "Operating Point No. 1" } ] }, "magnetic": {"coil": {"bobbin": "Basic", "functionalDescription":[{"name": "Primary", "numberTurns": 5, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.80 - Grade 1" }] }, "core": {"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } }, "manufacturerInfo": {"name": "", "reference": "Example" } } })");
    //     auto mas = OpenMagnetics::MasWrapper(masJson);

    //     auto coil = mas.get_mutable_magnetic().get_mutable_coil();
    //     auto core = mas.get_mutable_magnetic().get_mutable_core();
    //     auto inputs = mas.get_mutable_inputs();
    //     core.process_data();
    //     core.process_gap();
    //     coil.set_bobbin(OpenMagnetics::BobbinWrapper::create_quick_bobbin(core));
    //     coil.wind();

    //     OpenMagnetics::StrayCapacitance strayCapacitance;

    //     std::map<std::pair<std::string, std::string>, double> expectedValues = {
    //         {{"Primary parallel 0 turn 0", "Primary parallel 0 turn 1"}, 6e-12},
    //         {{"Primary parallel 0 turn 1", "Primary parallel 0 turn 2"}, 6e-12},
    //         {{"Primary parallel 0 turn 2", "Primary parallel 0 turn 3"}, 6e-12},
    //         {{"Primary parallel 0 turn 3", "Primary parallel 0 turn 4"}, 6e-12},
    //     };

    //     auto capacitanceAmongTurns = strayCapacitance.calculate_capacitance_among_turns(coil);
    //     CHECK(capacitanceAmongTurns.size() == 4);
    //     for (auto [keys, capacitance] : capacitanceAmongTurns) {
    //         // std::cout << "Capacitance between turn " << keys.first << " and turn " << keys.second << ": " << capacitance << std::endl;
    //         CHECK_CLOSE(expectedValues[keys], capacitance, expectedValues[keys] * maximumError);
    //     }
    // }

    // TEST(Test_Get_Capacitance_Among_Layers_Symmetrical_No_Interleaving) {
    //     settings->reset();
    //     json coilJson = json::parse(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 16, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 16, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] })");
    //     json coreJson = json::parse(R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })");

    //     auto core = OpenMagnetics::CoreWrapper(coreJson);
    //     auto coil = OpenMagnetics::CoilWrapper(coilJson);

    //     core.process_data();
    //     core.process_gap();
    //     coil.set_bobbin(OpenMagnetics::BobbinWrapper::create_quick_bobbin(core));
    //     coil.wind();

    //     OpenMagnetics::StrayCapacitance strayCapacitance;

    //     std::map<std::pair<std::string, std::string>, double> expectedValues = {
    //         {{"Primary section 0 layer 0", "Primary section 0 layer 0"}, 29e-12},
    //         {{"Primary section 0 layer 0", "Primary section 0 layer 1"}, 27e-12},
    //         {{"Primary section 0 layer 0", "Secondary section 0 layer 0"}, 0},
    //         {{"Primary section 0 layer 0", "Secondary section 0 layer 1"}, 0},
    //         {{"Primary section 0 layer 1", "Primary section 0 layer 1"}, 52e-12},
    //         {{"Primary section 0 layer 1", "Secondary section 0 layer 0"}, 25e-12},
    //         {{"Primary section 0 layer 1", "Secondary section 0 layer 1"}, 0},
    //         {{"Secondary section 0 layer 0", "Secondary section 0 layer 0"}, 59e-12},
    //         {{"Secondary section 0 layer 0", "Secondary section 0 layer 1"}, 33e-12},
    //         {{"Secondary section 0 layer 1", "Secondary section 0 layer 1"}, 35e-12},
    //     };

    //     auto capacitanceAmongLayers = strayCapacitance.calculate_capacitance_among_layers(coil);
    //     for (auto [keys, capacitance] : capacitanceAmongLayers) {
    //         std::cout << "Capacitance between " << keys.first << " and " << keys.second << ": " << capacitance << std::endl;
    //     }

    // }

    // TEST(Test_Get_Capacitance_Among_Windings_1_Winding_2_Turns_1_Parallel) {
    //     settings->reset();
    //     json coilJson = json::parse(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 2, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" } ] })");
    //     json coreJson = json::parse(R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })");

    //     auto core = OpenMagnetics::CoreWrapper(coreJson);
    //     auto coil = OpenMagnetics::CoilWrapper(coilJson);

    //     core.process_data();
    //     core.process_gap();
    //     coil.set_bobbin(OpenMagnetics::BobbinWrapper::create_quick_bobbin(core));
    //     coil.wind();

    //     OpenMagnetics::StrayCapacitance strayCapacitance;

    //     std::map<std::pair<std::string, std::string>, double> expectedValues = {
    //         {{"Primary", "Primary"}, 2.1e-12},
    //     };

    //     auto capacitanceAmongLayers = strayCapacitance.calculate_capacitance_among_windings(coil);
    //     for (auto [keys, capacitance] : capacitanceAmongLayers) {
    //         std::cout << "Capacitance between " << keys.first << " and " << keys.second << ": " << capacitance << std::endl;
    //     }

    // }

    // TEST(Test_Get_Capacitance_Among_Windings_1_Winding_3_Turns_1_Parallel) {
    //     settings->reset();
    //     json coilJson = json::parse(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 3, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" } ] })");
    //     json coreJson = json::parse(R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })");

    //     auto core = OpenMagnetics::CoreWrapper(coreJson);
    //     auto coil = OpenMagnetics::CoilWrapper(coilJson);

    //     core.process_data();
    //     core.process_gap();
    //     coil.set_bobbin(OpenMagnetics::BobbinWrapper::create_quick_bobbin(core));
    //     coil.wind();

    //     OpenMagnetics::StrayCapacitance strayCapacitance;

    //     std::map<std::pair<std::string, std::string>, double> expectedValues = {
    //         {{"Primary", "Primary"}, 1.9e-12},
    //     };

    //     auto capacitanceAmongLayers = strayCapacitance.calculate_capacitance_among_windings(coil);
    //     for (auto [keys, capacitance] : capacitanceAmongLayers) {
    //         std::cout << "Capacitance between " << keys.first << " and " << keys.second << ": " << capacitance << std::endl;
    //     }

    // }

    // TEST(Test_Get_Capacitance_Among_Windings_1_Winding_8_Turns_1_Parallel) {
    //     settings->reset();
    //     json coilJson = json::parse(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" } ] })");
    //     json coreJson = json::parse(R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })");

    //     auto core = OpenMagnetics::CoreWrapper(coreJson);
    //     auto coil = OpenMagnetics::CoilWrapper(coilJson);

    //     core.process_data();
    //     core.process_gap();
    //     coil.set_bobbin(OpenMagnetics::BobbinWrapper::create_quick_bobbin(core));
    //     coil.wind();

    //     OpenMagnetics::StrayCapacitance strayCapacitance;

    //     std::map<std::pair<std::string, std::string>, double> expectedValues = {
    //         {{"Primary", "Primary"}, 1.2e-12},
    //     };

    //     auto capacitanceAmongLayers = strayCapacitance.calculate_capacitance_among_windings(coil);
    //     for (auto [keys, capacitance] : capacitanceAmongLayers) {
    //         std::cout << "Capacitance between " << keys.first << " and " << keys.second << ": " << capacitance << std::endl;
    //     }

    // }

    // TEST(Test_Get_Capacitance_Among_Windings_1_Winding_16_Turns_1_Parallel_2_Layers) {
    //     settings->reset();
    //     json coilJson = json::parse(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 16, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" } ] })");
    //     json coreJson = json::parse(R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })");

    //     auto core = OpenMagnetics::CoreWrapper(coreJson);
    //     auto coil = OpenMagnetics::CoilWrapper(coilJson);

    //     core.process_data();
    //     core.process_gap();
    //     coil.set_bobbin(OpenMagnetics::BobbinWrapper::create_quick_bobbin(core));
    //     coil.wind();

    //     OpenMagnetics::StrayCapacitance strayCapacitance;

    //     std::map<std::pair<std::string, std::string>, double> expectedValues = {
    //         {{"Primary", "Primary"}, 14.9e-12},
    //     };

    //     auto capacitanceAmongLayers = strayCapacitance.calculate_capacitance_among_windings(coil);
    //     for (auto [keys, capacitance] : capacitanceAmongLayers) {
    //         std::cout << "Capacitance between " << keys.first << " and " << keys.second << ": " << capacitance << std::endl;
    //     }

    // }

    // TEST(Test_Get_Capacitance_Among_Windings_2_Windings_1_Turn_1_Parallels) {
    //     settings->reset();
    //     json coilJson = json::parse(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 1, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 1, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] })");
    //     json coreJson = json::parse(R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })");

    //     auto core = OpenMagnetics::CoreWrapper(coreJson);
    //     auto coil = OpenMagnetics::CoilWrapper(coilJson);

    //     core.process_data();
    //     core.process_gap();
    //     coil.set_bobbin(OpenMagnetics::BobbinWrapper::create_quick_bobbin(core));
    //     coil.wind();

    //     OpenMagnetics::StrayCapacitance strayCapacitance;

    //     std::map<std::pair<std::string, std::string>, double> expectedValues = {
    //         {{"Primary", "Primary"}, 1.5e-12},
    //         {{"Primary", "Secondary"}, 1.3e-12},
    //         {{"Secondary", "Secondary"}, 1.54e-12},
    //     };

    //     auto capacitanceAmongLayers = strayCapacitance.calculate_capacitance_among_windings(coil);
    //     for (auto [keys, capacitance] : capacitanceAmongLayers) {
    //         std::cout << "Capacitance between " << keys.first << " and " << keys.second << ": " << capacitance << std::endl;
    //     }

    // }

    // TEST(Test_Get_Capacitance_Among_Windings_2_Windings_2_Turns_1_Parallels) {
    //     settings->reset();
    //     json coilJson = json::parse(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 2, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 2, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] })");
    //     json coreJson = json::parse(R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })");

    //     auto core = OpenMagnetics::CoreWrapper(coreJson);
    //     auto coil = OpenMagnetics::CoilWrapper(coilJson);

    //     core.process_data();
    //     core.process_gap();
    //     coil.set_bobbin(OpenMagnetics::BobbinWrapper::create_quick_bobbin(core));
    //     coil.wind();

    //     OpenMagnetics::StrayCapacitance strayCapacitance;

    //     std::map<std::pair<std::string, std::string>, double> expectedValues = {
    //         {{"Primary", "Primary"}, 3.12e-12},
    //         {{"Primary", "Secondary"}, 2.2e-12},
    //         {{"Secondary", "Secondary"}, 3.3e-12},
    //     };

    //     // auto capacitanceAmongLayers = strayCapacitance.calculate_capacitance_among_turns(coil);
    //     auto capacitanceAmongLayers = strayCapacitance.calculate_capacitance_among_windings(coil);
    //     for (auto [keys, capacitance] : capacitanceAmongLayers) {
    //         std::cout << "Capacitance between " << keys.first << " and " << keys.second << ": " << capacitance << std::endl;
    //     }

    //     OpenMagnetics::MagneticWrapper magnetic;
    //     magnetic.set_core(core);
    //     magnetic.set_coil(coil);


    //     auto outFile = outputFilePath;
    //     outFile.append("Test_Get_Capacitance_Among_Windings_2_Windings_2_Turns_1_Parallels.svg");
    //     std::filesystem::remove(outFile);
    //     OpenMagnetics::Painter painter(outFile);
    //     painter.paint_core(magnetic);
    //     painter.paint_bobbin(magnetic);
    //     painter.paint_coil_turns(magnetic);
    //     painter.export_svg();

    // }


    // TEST(Test_Get_Capacitance_Among_Windings_2_Windings_8_Turns_1_Parallels) {
    //     settings->reset();
    //     json coilJson = json::parse(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] })");
    //     json coreJson = json::parse(R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })");

    //     auto core = OpenMagnetics::CoreWrapper(coreJson);
    //     auto coil = OpenMagnetics::CoilWrapper(coilJson);

    //     core.process_data();
    //     core.process_gap();
    //     coil.set_bobbin(OpenMagnetics::BobbinWrapper::create_quick_bobbin(core));
    //     coil.wind();

    //     OpenMagnetics::StrayCapacitance strayCapacitance;

    //     std::map<std::pair<std::string, std::string>, double> expectedValues = {
    //         {{"Primary", "Primary"}, 3.12e-12},
    //         {{"Primary", "Secondary"}, 2.2e-12},
    //         {{"Secondary", "Secondary"}, 3.3e-12},
    //     };

    //     // auto capacitanceAmongLayers = strayCapacitance.calculate_capacitance_among_turns(coil);
    //     auto capacitanceAmongLayers = strayCapacitance.calculate_capacitance_among_windings(coil);
    //     for (auto [keys, capacitance] : capacitanceAmongLayers) {
    //         std::cout << "Capacitance between " << keys.first << " and " << keys.second << ": " << capacitance << std::endl;
    //     }

    //     OpenMagnetics::MagneticWrapper magnetic;
    //     magnetic.set_core(core);
    //     magnetic.set_coil(coil);


    //     auto outFile = outputFilePath;
    //     outFile.append("Test_Get_Capacitance_Among_Windings_2_Windings_8_Turns_1_Parallels.svg");
    //     std::filesystem::remove(outFile);
    //     OpenMagnetics::Painter painter(outFile);
    //     painter.paint_core(magnetic);
    //     painter.paint_bobbin(magnetic);
    //     painter.paint_coil_turns(magnetic);
    //     painter.export_svg();
    // }

    TEST(Test_Get_Capacitance_Among_Layers_2_Windings_8_Turns_1_Parallels) {
        settings->reset();
        json coilJson = json::parse(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] })");
        json coreJson = json::parse(R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "RM 10/I", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })");

        auto core = OpenMagnetics::CoreWrapper(coreJson);
        auto coil = OpenMagnetics::CoilWrapper(coilJson);

        core.process_data();
        core.process_gap();
        coil.set_bobbin(OpenMagnetics::BobbinWrapper::create_quick_bobbin(core));
        coil.wind();

        OpenMagnetics::StrayCapacitance strayCapacitance;

        std::map<std::pair<std::string, std::string>, double> expectedValues = {
            {{"Primary", "Primary"}, 3.12e-12},
            {{"Primary", "Secondary"}, 2.2e-12},
            {{"Secondary", "Secondary"}, 3.3e-12},
        };

        // auto capacitanceAmongTurns = strayCapacitance.calculate_capacitance_among_turns(coil);
        // CHECK(capacitanceAmongTurns.size() == 4);
        // for (auto [keys, capacitance] : capacitanceAmongTurns) {
        //     std::cout << "Capacitance between turn " << keys.first << " and turn " << keys.second << ": " << capacitance << std::endl;
        // }

        auto capacitanceAmongLayers = strayCapacitance.calculate_capacitance_among_layers(coil);
        for (auto [keys, capacitance] : capacitanceAmongLayers) {
            std::cout << "Capacitance between " << keys.first << " and " << keys.second << ": " << capacitance << std::endl;
        }

        OpenMagnetics::MagneticWrapper magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);


        auto outFile = outputFilePath;
        outFile.append("Test_Get_Capacitance_Among_Windings_2_Windings_8_Turns_1_Parallels.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();

    }
    TEST(Test_Get_Capacitance_Among_Layers_1_Windings_6_Turns_1_Parallels) {
        settings->reset();
        json coilJson = json::parse(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 6, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" } ] })");
        json coreJson = json::parse(R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "PQ 32/12", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })");

        auto core = OpenMagnetics::CoreWrapper(coreJson);
        auto coil = OpenMagnetics::CoilWrapper(coilJson);

        core.process_data();
        core.process_gap();
        coil.set_bobbin(OpenMagnetics::BobbinWrapper::create_quick_bobbin(core));
        coil.wind();

        OpenMagnetics::StrayCapacitance strayCapacitance;

        std::map<std::pair<std::string, std::string>, double> expectedValues = {
            {{"Primary", "Primary"}, 3.12e-12},
            {{"Primary", "Secondary"}, 2.2e-12},
            {{"Secondary", "Secondary"}, 3.3e-12},
        };

        // auto capacitanceAmongTurns = strayCapacitance.calculate_capacitance_among_turns(coil);
        // CHECK(capacitanceAmongTurns.size() == 4);
        // for (auto [keys, capacitance] : capacitanceAmongTurns) {
        //     std::cout << "Capacitance between turn " << keys.first << " and turn " << keys.second << ": " << capacitance << std::endl;
        // }

        auto capacitanceAmongLayers = strayCapacitance.calculate_capacitance_among_layers(coil);
        for (auto [keys, capacitance] : capacitanceAmongLayers) {
            std::cout << "Capacitance between " << keys.first << " and " << keys.second << ": " << capacitance << std::endl;
        }

        OpenMagnetics::MagneticWrapper magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);


        auto outFile = outputFilePath;
        outFile.append("Test_Get_Capacitance_Among_Windings_2_Windings_8_Turns_1_Parallels.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();

    }


    // TEST(Test_Get_Surrounding_Turns_All) {
    //     settings->reset();
    //     std::string file_path = __FILE__;
    //     auto inventory_path = file_path.substr(0, file_path.rfind("/")).append("/testData/dont_worry_this_is_just_a_dummy_mas.json");
    //     std::ifstream jsonFile(inventory_path);
    //     json masJson;
    //     jsonFile >> masJson;

    //     auto mas = OpenMagnetics::MasWrapper(masJson);

    //     auto turns = mas.get_magnetic().get_coil().get_turns_description().value();

    //     size_t index = 0;
    //     for (auto turn : turns) {

    //         std::cout << turn.get_name() << std::endl;
    //         auto surroundingTurns = OpenMagnetics::StrayCapacitance::get_surrounding_turns(turn, turns);
    //         std::cout << surroundingTurns.size() << std::endl;
    //         std::vector<OpenMagnetics::Turn> newTurns = surroundingTurns;
    //         newTurns.push_back(turn);

    //         mas.get_mutable_magnetic().get_mutable_coil().set_turns_description(newTurns);

    //         auto outFile = outputFilePath;
    //         outFile.append("Test_Get_Surrounding_Turns_" + std::to_string(index) + ".svg");
    //         std::filesystem::remove(outFile);
    //         OpenMagnetics::Painter painter(outFile);
    //         painter.paint_core(mas.get_mutable_magnetic());
    //         painter.paint_bobbin(mas.get_mutable_magnetic());
    //         painter.paint_coil_turns(mas.get_mutable_magnetic());
    //         painter.export_svg();
    //         index++;
    //     }

    // }

}
