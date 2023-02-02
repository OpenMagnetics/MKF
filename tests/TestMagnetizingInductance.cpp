#include <UnitTest++.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <filesystem>
#include <nlohmann/json-schema.hpp>
#include "json.hpp"
#include "Utils.h"
#include "MagnetizingInductance.h"
#include <magic_enum.hpp>
#include <typeinfo>
#include "TestingUtils.h"


using json = nlohmann::json;


SUITE(MagnetizingInductance)
{
    double max_error = 0.05;
    void prepare_test_parameters(double dcCurrent,
                                 double ambientTemperature,
                                 double frequency,
                                 double numberTurns,
                                 double desiredMagnetizingInductance,
                                 std::vector<OpenMagnetics::CoreGap> gapping,
                                 std::string coreShape,
                                 std::string coreMaterial,
                                 OpenMagnetics::CoreWrapper& core,
                                 OpenMagnetics::WindingWrapper& winding,
                                 OpenMagnetics::InputsWrapper& inputs) {

        double dutyCycle = 0.5;
        double peakToPeak = 2;

        inputs = OpenMagnetics::InputsWrapper::create_quick_operation_point(frequency,
                                                                            desiredMagnetizingInductance,
                                                                            ambientTemperature,
                                                                            OpenMagnetics::WaveformLabel::SQUARE,
                                                                            peakToPeak,
                                                                            dutyCycle,
                                                                            dcCurrent);

        json primaryWindingJson = json();
        primaryWindingJson["isolationSide"] = OpenMagnetics::IsolationSide::PRIMARY;
        primaryWindingJson["name"] = "primary";
        primaryWindingJson["numberParallels"] = 1;
        primaryWindingJson["numberTurns"] = numberTurns;
        primaryWindingJson["wire"] = "Dummy";
        OpenMagnetics::WindingFunctionalDescription primaryWindingFunctionalDescription(primaryWindingJson);
        json windingFunctionalDescriptionJson = json::array();
        windingFunctionalDescriptionJson.push_back(primaryWindingJson);
        json windingJson = json();
        windingJson["functionalDescription"] = windingFunctionalDescriptionJson;
        OpenMagnetics::WindingWrapper windingAux(windingJson);
        winding = windingAux;

        core = OpenMagneticsTesting::get_core(coreShape, gapping, 1, coreMaterial);

    }

    TEST(Test_Inductance_Ferrite_Grinded)
    {
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 666;
        double frequency = 20000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.003);

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::WindingWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::map<std::string, std::string>({{"reluctance", "ZHANG"}}));

        double expectedValue = 23.3e-3;

        prepare_test_parameters(dcCurrent,
                                ambientTemperature,
                                frequency,
                                numberTurns,
                                -1,
                                gapping,
                                coreShape,
                                coreMaterial,
                                core,
                                winding,
                                inputs);

        double magnetizingInductance = magnetizing_inductance.get_inductance_from_number_turns_and_gapping(core, winding, inputs.get_operation_point(0));

        CHECK_CLOSE(expectedValue, magnetizingInductance, max_error * expectedValue);
    }


    TEST(Test_Inductance_Ferrite_Web)
    {
        // This tests checks that the operation is not crashing

        json coreData = json::parse("{\"functionalDescription\": {\"gapping\": [{\"area\": null, \"coordinates\": null, \"distanceClosestNormalSurface\": null, \"distanceClosestParallelSurface\": null, \"length\": 0.001, \"sectionDimensions\": null, \"shape\": null, \"type\": \"subtractive\"}, {\"area\": null, \"coordinates\": null, \"distanceClosestNormalSurface\": null, \"distanceClosestParallelSurface\": null, \"length\": 1e-05, \"sectionDimensions\": null, \"shape\": null, \"type\": \"residual\"}, {\"area\": null, \"coordinates\": null, \"distanceClosestNormalSurface\": null, \"distanceClosestParallelSurface\": null, \"length\": 1e-05, \"sectionDimensions\": null, \"shape\": null, \"type\": \"residual\"}], \"material\": \"3C97\", \"name\": \"My Core\", \"numberStacks\": 1, \"shape\": {\"aliases\": [], \"dimensions\": {\"A\": 0.0391, \"B\": 0.0198, \"C\": 0.0125, \"D\": 0.0146, \"E\": 0.030100000000000002, \"F\": 0.0125, \"G\": 0.0, \"H\": 0.0}, \"family\": \"etd\", \"familySubtype\": \"1\", \"magneticCircuit\": null, \"name\": \"ETD 39/20/13\", \"type\": \"standard\"}, \"type\": \"two-piece set\"}}");
        json windingData = json::parse("{\"functionalDescription\": [{\"isolationSide\": \"primary\", \"name\": \"Primary\", \"numberParallels\": 1, \"numberTurns\": 1, \"wire\": \"Dummy\"}], \"layersDescription\": null, \"sectionsDescription\": null, \"turnsDescription\": null}");
        json operationPointData = json::parse("{\"conditions\": {\"ambientRelativeHumidity\": null, \"ambientTemperature\": 25.0, \"cooling\": null, \"name\": null}, \"excitationsPerWinding\": [{\"current\": {\"harmonics\": null, \"processed\": null, \"waveform\": {\"ancillaryLabel\": null, \"data\": [-5.0, 5.0, -5.0], \"numberPeriods\": null, \"time\": [0.0, 2.5e-06, 1e-05]}}, \"frequency\": 100000.0, \"magneticField\": null, \"magneticFluxDensity\": null, \"magnetizingCurrent\": null, \"name\": \"My Operation Point\", \"voltage\": {\"harmonics\": null, \"processed\": null, \"waveform\": {\"ancillaryLabel\": null, \"data\": [7.5, 7.5, -2.5, -2.5, 7.5], \"numberPeriods\": null, \"time\": [0.0, 2.5e-06, 2.5e-06, 1e-05, 1e-05]}}}],\"name\": null}");

        OpenMagnetics::CoreWrapper core(coreData);
        OpenMagnetics::WindingWrapper winding(windingData);
        OpenMagnetics::OperationPoint operationPoint(operationPointData);
        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::map<std::string, std::string>({{"reluctance", "ZHANG"}}));
        magnetizing_inductance.get_inductance_from_number_turns_and_gapping(core, winding, operationPoint);
    }

    TEST(Test_Inductance_Ferrite_Spacer)
    {
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 666;
        double frequency = 20000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_spacer_gap(0.003);
        double expectedValue = 13.5e-3;

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::WindingWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::map<std::string, std::string>({{"reluctance", "ZHANG"}}));

        prepare_test_parameters(dcCurrent,
                                ambientTemperature,
                                frequency,
                                numberTurns,
                                -1,
                                gapping,
                                coreShape,
                                coreMaterial,
                                core,
                                winding,
                                inputs);

        double magnetizingInductance = magnetizing_inductance.get_inductance_from_number_turns_and_gapping(core, winding, inputs.get_operation_point(0));

        CHECK_CLOSE(expectedValue, magnetizingInductance, max_error * expectedValue);
    }

    TEST(Test_Inductance_Ferrite_Grinded_Few_Turns)
    {
        double dcCurrent = 0;
        double ambientTemperature = 42;
        double numberTurns = 9;
        double frequency = 100000;
        std::string coreShape = "E 47/20/16";
        std::string coreMaterial = "N87";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0004);
        double expectedValue = 63e-6;

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::WindingWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::map<std::string, std::string>({{"reluctance", "ZHANG"}}));

        prepare_test_parameters(dcCurrent,
                                ambientTemperature,
                                frequency,
                                numberTurns,
                                -1,
                                gapping,
                                coreShape,
                                coreMaterial,
                                core,
                                winding,
                                inputs);

        double magnetizingInductance = magnetizing_inductance.get_inductance_from_number_turns_and_gapping(core, winding, inputs.get_operation_point(0));

        CHECK_CLOSE(expectedValue, magnetizingInductance, max_error * expectedValue);
    }

    TEST(Test_Inductance_Powder)
    {
        double dcCurrent = 96;
        double ambientTemperature = 25;
        double numberTurns = 13;
        double frequency = 68000;
        std::string coreShape = "E 42/21/15";
        std::string coreMaterial = "Edge 60";
        auto gapping = OpenMagneticsTesting::get_residual_gap();
        double expectedValue = 15.7e-6;

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::WindingWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::map<std::string, std::string>({{"reluctance", "ZHANG"}}));

        prepare_test_parameters(dcCurrent,
                                ambientTemperature,
                                frequency,
                                numberTurns,
                                -1,
                                gapping,
                                coreShape,
                                coreMaterial,
                                core,
                                winding,
                                inputs);

        double magnetizingInductance = magnetizing_inductance.get_inductance_from_number_turns_and_gapping(core, winding, inputs.get_operation_point(0)) ;

        CHECK_CLOSE(expectedValue, magnetizingInductance, max_error * expectedValue);
    }

    TEST(Test_NumberTurns_Ferrite_Grinded)
    {
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 23.3e-3;
        double frequency = 20000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.003);

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::WindingWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::map<std::string, std::string>({{"reluctance", "ZHANG"}}));

        double expectedValue = 666;

        prepare_test_parameters(dcCurrent,
                                ambientTemperature,
                                frequency,
                                -1,
                                desiredMagnetizingInductance,
                                gapping,
                                coreShape,
                                coreMaterial,
                                core,
                                winding,
                                inputs);

        double numberTurns = magnetizing_inductance.get_number_turns_from_gapping_and_inductance(core, inputs);

        CHECK_CLOSE(expectedValue, numberTurns, max_error * expectedValue);
    }

    TEST(Test_NumberTurns_Powder)
    {
        double dcCurrent = 96;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 15.7e-6;
        double frequency = 68000;
        std::string coreShape = "E 42/21/15";
        std::string coreMaterial = "Edge 60";
        auto gapping = OpenMagneticsTesting::get_residual_gap();

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::WindingWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::map<std::string, std::string>({{"reluctance", "ZHANG"}}));

        double expectedValue = 13;

        prepare_test_parameters(dcCurrent,
                                ambientTemperature,
                                frequency,
                                -1,
                                desiredMagnetizingInductance,
                                gapping,
                                coreShape,
                                coreMaterial,
                                core,
                                winding,
                                inputs);

        double numberTurns = magnetizing_inductance.get_number_turns_from_gapping_and_inductance(core, inputs);

        CHECK_CLOSE(expectedValue, numberTurns, max_error * expectedValue);
    }

    TEST(Test_Gapping_Ferrite_Grinded)
    {
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 23.3e-3;
        double numberTurns = 666;
        double frequency = 20000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::WindingWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::map<std::string, std::string>({{"reluctance", "ZHANG"}}));

        double expectedValue = 0.003;

        prepare_test_parameters(dcCurrent,
                                ambientTemperature,
                                frequency,
                                numberTurns,
                                desiredMagnetizingInductance,
                                {},
                                coreShape,
                                coreMaterial,
                                core,
                                winding,
                                inputs);

        auto gapping = magnetizing_inductance.get_gapping_from_number_turns_and_inductance(core, winding, inputs, OpenMagnetics::GappingType::GRINDED);

        CHECK_CLOSE(expectedValue, gapping[0].get_length(), max_error * expectedValue);
    }

    TEST(Test_Gapping_U_Shape_Ferrite_Grinded)
    {
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 23.3e-3;
        double numberTurns = 666;
        double frequency = 20000;
        std::string coreShape = "U 26/22/16";
        std::string coreMaterial = "3C97";

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::WindingWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::map<std::string, std::string>({{"reluctance", "ZHANG"}}));

        double expectedValue = 0.0066;

        prepare_test_parameters(dcCurrent,
                                ambientTemperature,
                                frequency,
                                numberTurns,
                                desiredMagnetizingInductance,
                                {},
                                coreShape,
                                coreMaterial,
                                core,
                                winding,
                                inputs);

        auto gapping = magnetizing_inductance.get_gapping_from_number_turns_and_inductance(core, winding, inputs, OpenMagnetics::GappingType::GRINDED);

        CHECK_CLOSE(expectedValue, gapping[0].get_length(), max_error * expectedValue);
    }

    TEST(Test_Gapping_Ferrite_Distributed)
    {
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 23.3e-3;
        double numberTurns = 666;
        double frequency = 20000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::WindingWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::map<std::string, std::string>({{"reluctance", "ZHANG"}}));

        double expectedValue = 0.0004;

        prepare_test_parameters(dcCurrent,
                                ambientTemperature,
                                frequency,
                                numberTurns,
                                desiredMagnetizingInductance,
                                {},
                                coreShape,
                                coreMaterial,
                                core,
                                winding,
                                inputs);

        auto gapping = magnetizing_inductance.get_gapping_from_number_turns_and_inductance(core, winding, inputs, OpenMagnetics::GappingType::DISTRIBUTED);

        CHECK_CLOSE(expectedValue, gapping[0].get_length(), max_error * expectedValue);
        CHECK_EQUAL(7UL, gapping.size());
    }

    TEST(Test_Gapping_Ferrite_Distributed_More_Gap_Precision)
    {
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 23.3e-3;
        double numberTurns = 666;
        double frequency = 20000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::WindingWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::map<std::string, std::string>({{"reluctance", "ZHANG"}}));

        double expectedValue = 0.00039;

        prepare_test_parameters(dcCurrent,
                                ambientTemperature,
                                frequency,
                                numberTurns,
                                desiredMagnetizingInductance,
                                {},
                                coreShape,
                                coreMaterial,
                                core,
                                winding,
                                inputs);

        auto gapping = magnetizing_inductance.get_gapping_from_number_turns_and_inductance(core, winding, inputs, OpenMagnetics::GappingType::DISTRIBUTED, 5);

        CHECK_EQUAL(expectedValue, gapping[0].get_length());
        CHECK_EQUAL(7UL, gapping.size());
    }
}