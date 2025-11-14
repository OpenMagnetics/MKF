#include "support/Painter.h"
#include "physical_models/StrayCapacitance.h"
#include "constructive_models/Mas.h"
#include "support/Utils.h"
#include "processors/Inputs.h"
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


SUITE(StrayCapacitance) {
    double maximumError = 0.4;
    auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
    bool plot = false;

    TEST(Test_Get_Capacitance_Among_Windings_2_Windings_1_Turn_1_Parallels) {
        settings->reset();
        json coilJson = json::parse(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 1, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 1, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] })");
        json coreJson = json::parse(R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })");

        auto core = OpenMagnetics::Core(coreJson);
        auto coil = OpenMagnetics::Coil(coilJson);

        core.process_data();
        core.process_gap();
        coil.set_bobbin(OpenMagnetics::Bobbin::create_quick_bobbin(core));
        coil.wind();

        StrayCapacitance strayCapacitance;

        std::map<std::pair<std::string, std::string>, double> expectedValues = {
            {{"Primary", "Primary"}, 2.8708e-12},
            {{"Primary", "Secondary"}, -2.4259e-12},
            {{"Secondary", "Primary"}, -2.4259e-12},
            {{"Secondary", "Secondary"}, 2.887e-12},
        };

        auto maxwellCapacitanceMatrix = strayCapacitance.calculate_maxwell_capacitance_matrix(coil);
        CHECK(maxwellCapacitanceMatrix.size() == 4);
        for (auto [keys, capacitance] : maxwellCapacitanceMatrix) {
            CHECK_CLOSE(expectedValues[keys], capacitance, fabs(expectedValues[keys]) * maximumError);
        }
    }

    TEST(Test_Get_Capacitance_Among_Windings_1_Windings_8_Turns_1_Parallels) {
        settings->reset();
        json coilJson = json::parse(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" } ] })");
        json coreJson = json::parse(R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "RM 10/I", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })");

        auto core = OpenMagnetics::Core(coreJson);
        auto coil = OpenMagnetics::Coil(coilJson);

        core.process_data();
        core.process_gap();
        coil.set_bobbin(OpenMagnetics::Bobbin::create_quick_bobbin(core));
        coil.wind();

        StrayCapacitance strayCapacitance;

        std::map<std::pair<std::string, std::string>, double> expectedValues = {
            {{"Primary", "Primary"}, 10.166e-12 - 9.8337e-12},
        };


        auto maxwellCapacitanceMatrix = strayCapacitance.calculate_maxwell_capacitance_matrix(coil);
        CHECK(maxwellCapacitanceMatrix.size() == 1);
        for (auto [keys, capacitance] : maxwellCapacitanceMatrix) {

            CHECK_CLOSE(expectedValues[keys], capacitance, fabs(expectedValues[keys]) * maximumError);
        }

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);


        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Get_Capacitance_Among_Windings_1_Windings_8_Turns_1_Parallels.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
    }

    TEST(Test_Get_Capacitance_Among_Windings_2_Windings_8_Turns_1_Parallels) {
        settings->reset();
        json coilJson = json::parse(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] })");
        json coreJson = json::parse(R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "RM 10/I", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })");

        auto core = OpenMagnetics::Core(coreJson);
        auto coil = OpenMagnetics::Coil(coilJson);

        core.process_data();
        core.process_gap();
        coil.set_bobbin(OpenMagnetics::Bobbin::create_quick_bobbin(core));
        coil.wind();

        StrayCapacitance strayCapacitance;

        std::map<std::pair<std::string, std::string>, double> expectedValues = {
            {{"Primary", "Primary"}, 10.166e-12},
            {{"Primary", "Secondary"}, -9.8337e-12},
            {{"Secondary", "Primary"}, -9.8337e-12},
            {{"Secondary", "Secondary"}, 10.199e-12},
        };


        auto maxwellCapacitanceMatrix = strayCapacitance.calculate_maxwell_capacitance_matrix(coil);
        CHECK(maxwellCapacitanceMatrix.size() == 4);
        for (auto [keys, capacitance] : maxwellCapacitanceMatrix) {
            CHECK_CLOSE(expectedValues[keys], capacitance, fabs(expectedValues[keys]) * maximumError);
        }

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);


        auto outFile = outputFilePath;
        outFile.append("Test_Get_Capacitance_Among_Windings_2_Windings_8_Turns_1_Parallels.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();

    }
    TEST(Test_Get_Capacitance_Among_Windings_2_Windings_16_8_Turns_1_Parallels) {
        settings->reset();
        json coilJson = json::parse(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 16, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] })");
        json coreJson = json::parse(R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "RM 10/I", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })");

        auto core = OpenMagnetics::Core(coreJson);
        auto coil = OpenMagnetics::Coil(coilJson);

        core.process_data();
        core.process_gap();
        coil.set_bobbin(OpenMagnetics::Bobbin::create_quick_bobbin(core));
        coil.wind();

        StrayCapacitance strayCapacitance;

        std::map<std::pair<std::string, std::string>, double> expectedValues = {
            {{"Primary", "Primary"}, 14.777e-12},
            {{"Primary", "Secondary"}, -6.6715e-12},
            {{"Secondary", "Primary"}, -6.6715e-12},
            {{"Secondary", "Secondary"}, 12.791e-12},
        };

        auto maxwellCapacitanceMatrix = strayCapacitance.calculate_maxwell_capacitance_matrix(coil);
        CHECK(maxwellCapacitanceMatrix.size() == 4);
        for (auto [keys, capacitance] : maxwellCapacitanceMatrix) {
            CHECK_CLOSE(expectedValues[keys], capacitance, fabs(expectedValues[keys]) * maximumError);
        }

    }

    TEST(Test_Get_Capacitance_Among_Windings_3_Windings_8_Turns_1_Parallels) {
        settings->reset();
        json coilJson = json::parse(R"({"bobbin": "Dummy", "functionalDescription":[{"name": "Primary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" }, {"name": "Tertiary", "numberTurns": 8, "numberParallels": 1, "isolationSide": "tertiary", "wire": "Round 1.00 - Grade 1" } ] })");
        json coreJson = json::parse(R"({"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "RM 10/I", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } })");

        auto core = OpenMagnetics::Core(coreJson);
        auto coil = OpenMagnetics::Coil(coilJson);

        core.process_data();
        core.process_gap();
        coil.set_bobbin(OpenMagnetics::Bobbin::create_quick_bobbin(core));
        coil.wind();

        StrayCapacitance strayCapacitance;

        std::map<std::pair<std::string, std::string>, double> expectedValues = {
            {{"Primary", "Primary"}, 15.437e-12},
            {{"Primary", "Secondary"}, -15.115e-12},
            {{"Primary", "Tertiary"}, -0.12966e-12},
            {{"Secondary", "Primary"}, -15.115e-12},
            {{"Secondary", "Secondary"}, 27.457e-12},
            {{"Secondary", "Tertiary"}, -12.21e-12},
            {{"Tertiary", "Primary"}, -0.12966e-12},
            {{"Tertiary", "Secondary"}, -12.21e-12},
            {{"Tertiary", "Tertiary"}, 12.791e-12},
        };

        auto maxwellCapacitanceMatrix = strayCapacitance.calculate_maxwell_capacitance_matrix(coil);
        CHECK(maxwellCapacitanceMatrix.size() == 9);
        for (auto [keys, capacitance] : maxwellCapacitanceMatrix) {
            if (capacitance != 0) {
                CHECK_CLOSE(expectedValues[keys], capacitance, fabs(expectedValues[keys]) * maximumError);
            }
            else {
                CHECK_CLOSE(expectedValues[keys], capacitance, 1e-12);
            }
        }

    }

    TEST(Test_One_Versus_Two_Windings) {
        {
            std::vector<int64_t> numberTurns = {110};
            std::vector<int64_t> numberParallels = {1};
            std::string shapeName = "T 12.5/7.5/5";
            std::vector<OpenMagnetics::Wire> wires;
            auto wire = find_wire_by_name("Round 0.15 - Grade 1");
            wires.push_back(wire);
            wires.push_back(wire);

            WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
            WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
            CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
            CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                             numberParallels,
                                                             shapeName,
                                                             1,
                                                             windingOrientation,
                                                             layersOrientation,
                                                             turnsAlignment,
                                                             sectionsAlignment,
                                                             wires,
                                                             false);

            int64_t numberStacks = 1;
            std::string coreMaterial = "A07";
            std::vector<CoreGap> gapping = {};
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            StrayCapacitance strayCapacitance;

            std::map<std::pair<std::string, std::string>, double> expectedValues = {
                {{"winding 0", "winding 0"}, 3.15617e-12 -3.14887e-12},
            };

            auto maxwellCapacitanceMatrix = strayCapacitance.calculate_maxwell_capacitance_matrix(coil);
            CHECK(maxwellCapacitanceMatrix.size() == 1);
            for (auto [keys, capacitance] : maxwellCapacitanceMatrix) {
                CHECK_CLOSE(expectedValues[keys], capacitance, fabs(expectedValues[keys]) * maximumError);
            }


            if (plot) {
                auto outFile = outputFilePath;
                outFile.append("Test_One_Versus_Two_Windings_One.svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile);
                painter.paint_core(magnetic);
                painter.paint_bobbin(magnetic);
                painter.paint_coil_turns(magnetic);
                painter.export_svg();
            }
        }
        {
            std::vector<int64_t> numberTurns = {110, 110};
            std::vector<int64_t> numberParallels = {1, 1};
            std::string shapeName = "T 12.5/7.5/5";
            std::vector<OpenMagnetics::Wire> wires;
            auto wire = find_wire_by_name("Round 0.15 - Grade 1");
            wires.push_back(wire);
            wires.push_back(wire);

            WindingOrientation windingOrientation = WindingOrientation::OVERLAPPING;
            WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
            CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
            CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                             numberParallels,
                                                             shapeName,
                                                             1,
                                                             windingOrientation,
                                                             layersOrientation,
                                                             turnsAlignment,
                                                             sectionsAlignment,
                                                             wires,
                                                             false);

            int64_t numberStacks = 1;
            std::string coreMaterial = "A07";
            std::vector<CoreGap> gapping = {};
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            StrayCapacitance strayCapacitance;

            std::map<std::pair<std::string, std::string>, double> expectedValues = {
                {{"winding 0", "winding 0"}, 3.15617e-12},
                {{"winding 0", "winding 1"}, -3.14887e-12},
                {{"winding 1", "winding 0"}, -3.14887e-12},
                {{"winding 1", "winding 1"}, 3.15674e-12},
            };

            auto maxwellCapacitanceMatrix = strayCapacitance.calculate_maxwell_capacitance_matrix(coil);
            CHECK(maxwellCapacitanceMatrix.size() == 4);
            for (auto [keys, capacitance] : maxwellCapacitanceMatrix) {

                CHECK_CLOSE(expectedValues[keys], capacitance, fabs(expectedValues[keys]) * maximumError);
            }

            if (plot) {
                auto outFile = outputFilePath;
                outFile.append("Test_One_Versus_Two_Windings_Two.svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile);
                painter.paint_core(magnetic);
                painter.paint_bobbin(magnetic);
                painter.paint_coil_turns(magnetic);
                painter.export_svg();
            }
        }

    }

    // TEST(Test_Get_Capacitance_Large) {
    //     settings->reset();
    //     json coilJson = json::parse(R"({"_interleavingLevel": 1, "_windingOrientation": "overlapping", "_layersOrientation": "overlapping", "_turnsAlignment": "spread", "_sectionAlignment": "inner or top", "bobbin": {"processedDescription": {"columnDepth": 0.0062, "columnShape": "rectangular", "columnThickness": 0.0011999999999999997, "columnWidth": 0.0062, "coordinates": [0, 0, 0 ], "wallThickness": 0.0011999999999999997, "windingWindows": [{"coordinates": [0.00935, 0, 0 ], "height": 0.022600000000000002, "width": 0.006300000000000001 } ] } }, "functionalDescription": [{"connections": [{"pinName": "4", "type": "Pin" }, {"pinName": "5", "type": "Pin" }, {"pinName": "6", "type": "Pin" }, {"pinName": "7", "type": "Pin" }, {"pinName": "8", "type": "Pin" }, {"pinName": "9", "type": "Pin" }, {"pinName": "10", "type": "Pin" }, {"pinName": "11", "type": "Pin" } ], "isolationSide": "primary", "name": "primary", "numberParallels": 16, "numberTurns": 4, "wire": "Round 25.0 - Single Build" }, {"connections": [{"pinName": "13", "type": "Pin" }, {"pinName": "3", "type": "Pin" } ], "isolationSide": "secondary", "name": "secondary", "numberParallels": 2, "numberTurns": 40, "wire": {"conductingDiameter": {"maximum": 0.000457, "minimum": 0.00045000000000000004, "nominal": 0.000455 }, "material": "copper", "outerDiameter": {"maximum": 0.000516, "minimum": 0.000495, "nominal": 0.000505 }, "coating": {"breakdownVoltage": 2370, "grade": 2, "type": "enamelled" }, "manufacturerInfo": {"name": "Elektrisola" }, "name": "Round 25.0 - Heavy Build", "numberConductors": 1, "standard": "NEMA MW 1000 C", "standardName": "25.0", "type": "round" } } ] })");
    //     json coreJson = json::parse(R"({"functionalDescription": {"gapping": [{"area": 0.0001, "coordinates": [0, 0.00055, 0 ], "distanceClosestNormalSurface": 0.01195, "distanceClosestParallelSurface": 0.007500000000000001, "length": 0.0011, "sectionDimensions": [0.01, 0.01 ], "shape": "rectangular", "type": "subtractive" }, {"area": 0.000051, "coordinates": [0.015001, 0, 0 ], "distanceClosestNormalSurface": 0.0125, "distanceClosestParallelSurface": 0.007500000000000001, "length": 0.000005, "sectionDimensions": [0.005001, 0.01 ], "shape": "rectangular", "type": "residual" }, {"area": 0.000051, "coordinates": [-0.015001, 0, 0 ], "distanceClosestNormalSurface": 0.0125, "distanceClosestParallelSurface": 0.007500000000000001, "length": 0.000005, "sectionDimensions": [0.005001, 0.01 ], "shape": "rectangular", "type": "residual" } ], "material": "TP4A", "numberStacks": 1, "shape": {"aliases": ["E 35" ], "dimensions": {"A": {"maximum": 0.0355, "minimum": 0.0345 }, "B": {"maximum": 0.01775, "minimum": 0.01725 }, "C": {"maximum": 0.0103, "minimum": 0.0097 }, "D": {"maximum": 0.01275, "minimum": 0.01225 }, "E": {"maximum": 0.0255, "minimum": 0.0245 }, "F": {"maximum": 0.0103, "minimum": 0.0097 } }, "family": "e", "magneticCircuit": "open", "name": "E 35/18/10", "type": "standard" }, "type": "two-piece set" }, "processedDescription": {"columns": [{"area": 0.0001, "coordinates": [0, 0, 0 ], "depth": 0.01, "height": 0.025, "shape": "rectangular", "type": "central", "width": 0.01 }, {"area": 0.000051, "coordinates": [0.015001, 0, 0 ], "depth": 0.01, "height": 0.025, "shape": "rectangular", "type": "lateral", "width": 0.005001 }, {"area": 0.000051, "coordinates": [-0.015001, 0, 0 ], "depth": 0.01, "height": 0.025, "shape": "rectangular", "type": "lateral", "width": 0.005001 } ], "depth": 0.01, "effectiveParameters": {"effectiveArea": 0.00010000000000000002, "effectiveLength": 0.08070796326794898, "effectiveVolume": 0.000008070796326794898, "minimumArea": 0.0001 }, "height": 0.035, "width": 0.035, "windingWindows": [{"area": 0.00018750000000000003, "coordinates": [0.005, 0 ], "height": 0.025, "width": 0.007500000000000001 } ] } })");

    //     auto core = OpenMagnetics::Core(coreJson);
    //     auto coil = OpenMagnetics::Coil(coilJson);

    //     core.process_data();
    //     core.process_gap();
    //     coil.set_bobbin(OpenMagnetics::Bobbin::create_quick_bobbin(core));
    //     coil.wind();

    //     StrayCapacitance strayCapacitance;

    //     std::map<std::pair<std::string, std::string>, double> expectedValues = {
    //         {{"Primary", "Primary"}, 3.12e-12},
    //         {{"Primary", "Secondary"}, 2.2e-12},
    //         {{"Secondary", "Secondary"}, 3.3e-12},
    //     };

    //     // auto capacitanceAmongTurns = strayCapacitance.calculate_capacitance_among_turns(coil);
    //     // CHECK(capacitanceAmongTurns.size() == 4);
    //     // for (auto [keys, capacitance] : capacitanceAmongTurns) {
    //     //     std::cout << "Capacitance between turn " << keys.first << " and turn " << keys.second << ": " << capacitance << std::endl;
    //     // }

    //     auto capacitanceAmongLayers = strayCapacitance.calculate_capacitance_among_windings(coil);
    //     for (auto [keys, capacitance] : capacitanceAmongLayers) {
    //         std::cout << "Capacitance between " << keys.first << " and " << keys.second << ": " << capacitance << std::endl;
    //     }

    //     OpenMagnetics::Magnetic magnetic;
    //     magnetic.set_core(core);
    //     magnetic.set_coil(coil);


    //     auto outFile = outputFilePath;
    //     outFile.append("Test_Get_Capacitance_Among_Windings_2_Windings_8_Turns_1_Parallels.svg");
    //     std::filesystem::remove(outFile);
    //     Painter painter(outFile);
    //     painter.paint_core(magnetic);
    //     painter.paint_bobbin(magnetic);
    //     painter.paint_coil_turns(magnetic);
    //     painter.export_svg();

    // }


    TEST(Test_Get_Surrounding_Turns_All) {
        settings->reset();
        json coilJson = json::parse(R"({"_interleavingLevel": 1, "_windingOrientation": "overlapping", "_layersOrientation": "overlapping", "_turnsAlignment": "spread", "_sectionAlignment": "inner or top", "bobbin": {"processedDescription": {"columnDepth": 0.0062, "columnShape": "rectangular", "columnThickness": 0.0011999999999999997, "columnWidth": 0.0062, "coordinates": [0, 0, 0 ], "wallThickness": 0.0011999999999999997, "windingWindows": [{"coordinates": [0.00935, 0, 0 ], "height": 0.022600000000000002, "width": 0.006300000000000001 } ] } }, "functionalDescription": [{"connections": [{"pinName": "4", "type": "Pin" }, {"pinName": "5", "type": "Pin" }, {"pinName": "6", "type": "Pin" }, {"pinName": "7", "type": "Pin" }, {"pinName": "8", "type": "Pin" }, {"pinName": "9", "type": "Pin" }, {"pinName": "10", "type": "Pin" }, {"pinName": "11", "type": "Pin" } ], "isolationSide": "primary", "name": "primary", "numberParallels": 2, "numberTurns": 16, "wire": "Round 19.0 - Heavy Build" }, {"connections": [{"pinName": "13", "type": "Pin" }, {"pinName": "3", "type": "Pin" } ], "isolationSide": "secondary", "name": "secondary", "numberParallels": 2, "numberTurns": 40, "wire": {"conductingDiameter": {"maximum": 0.000457, "minimum": 0.00045000000000000004, "nominal": 0.000455 }, "material": "copper", "outerDiameter": {"maximum": 0.000516, "minimum": 0.000495, "nominal": 0.000505 }, "coating": {"breakdownVoltage": 2370, "grade": 2, "type": "enamelled" }, "manufacturerInfo": {"name": "Elektrisola" }, "name": "Round 25.0 - Heavy Build", "numberConductors": 1, "standard": "NEMA MW 1000 C", "standardName": "25.0", "type": "round" } } ] })");
        json coreJson = json::parse(R"({"functionalDescription": {"gapping": [{"area": 0.0001, "coordinates": [0, 0.00055, 0 ], "distanceClosestNormalSurface": 0.01195, "distanceClosestParallelSurface": 0.007500000000000001, "length": 0.0011, "sectionDimensions": [0.01, 0.01 ], "shape": "rectangular", "type": "subtractive" }, {"area": 0.000051, "coordinates": [0.015001, 0, 0 ], "distanceClosestNormalSurface": 0.0125, "distanceClosestParallelSurface": 0.007500000000000001, "length": 0.000005, "sectionDimensions": [0.005001, 0.01 ], "shape": "rectangular", "type": "residual" }, {"area": 0.000051, "coordinates": [-0.015001, 0, 0 ], "distanceClosestNormalSurface": 0.0125, "distanceClosestParallelSurface": 0.007500000000000001, "length": 0.000005, "sectionDimensions": [0.005001, 0.01 ], "shape": "rectangular", "type": "residual" } ], "material": "TP4A", "numberStacks": 1, "shape": {"aliases": ["E 35" ], "dimensions": {"A": {"maximum": 0.0355, "minimum": 0.0345 }, "B": {"maximum": 0.01775, "minimum": 0.01725 }, "C": {"maximum": 0.0103, "minimum": 0.0097 }, "D": {"maximum": 0.01275, "minimum": 0.01225 }, "E": {"maximum": 0.0255, "minimum": 0.0245 }, "F": {"maximum": 0.0103, "minimum": 0.0097 } }, "family": "e", "magneticCircuit": "open", "name": "E 35/18/10", "type": "standard" }, "type": "two-piece set" }, "processedDescription": {"columns": [{"area": 0.0001, "coordinates": [0, 0, 0 ], "depth": 0.01, "height": 0.025, "shape": "rectangular", "type": "central", "width": 0.01 }, {"area": 0.000051, "coordinates": [0.015001, 0, 0 ], "depth": 0.01, "height": 0.025, "shape": "rectangular", "type": "lateral", "width": 0.005001 }, {"area": 0.000051, "coordinates": [-0.015001, 0, 0 ], "depth": 0.01, "height": 0.025, "shape": "rectangular", "type": "lateral", "width": 0.005001 } ], "depth": 0.01, "effectiveParameters": {"effectiveArea": 0.00010000000000000002, "effectiveLength": 0.08070796326794898, "effectiveVolume": 0.000008070796326794898, "minimumArea": 0.0001 }, "height": 0.035, "width": 0.035, "windingWindows": [{"area": 0.00018750000000000003, "coordinates": [0.005, 0 ], "height": 0.025, "width": 0.007500000000000001 } ] } })");

        auto core = OpenMagnetics::Core(coreJson);
        auto coil = OpenMagnetics::Coil(coilJson);

        core.process_data();
        core.process_gap();
        coil.set_bobbin(OpenMagnetics::Bobbin::create_quick_bobbin(core));
        coil.wind();

        auto turns = coil.get_turns_description().value();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        size_t index = 0;
        for (auto turn : turns) {

            auto surroundingTurns = StrayCapacitance::get_surrounding_turns(turn, turns);
            std::vector<Turn> newTurns;
            for (auto [surroundingTurn, surroundingTurnIndex] : surroundingTurns) {
                newTurns.push_back(surroundingTurn);
            }
            newTurns.push_back(turn);

            // magnetic.get_mutable_coil().set_turns_description(newTurns);

            auto outFile = outputFilePath;
            outFile.append("Test_Get_Surrounding_Turns_" + std::to_string(index) + ".svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            index++;
            break;
        }

    }

}
