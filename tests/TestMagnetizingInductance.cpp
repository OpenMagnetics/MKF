#include <source_location>
#include "processors/CircuitSimulatorInterface.h"
#include "physical_models/MagnetizingInductance.h"
#include "constructive_models/Bobbin.h"
#include "TestingUtils.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <typeinfo>
#include <vector>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace MAS;
using namespace OpenMagnetics;

using json = nlohmann::json;

namespace { 
    double max_error = 0.05;
    void prepare_test_parameters(double dcCurrent, double ambientTemperature, double frequency, double numberTurns,
                                 double desiredMagnetizingInductance, std::vector<CoreGap> gapping,
                                 std::string coreShape, std::string coreMaterial, Core& core,
                                 OpenMagnetics::Coil& winding,  OpenMagnetics::Inputs& inputs,
                                 double peakToPeak = 20, int numberStacks = 1) {
        double dutyCycle = 0.5;

        inputs = OpenMagnetics::Inputs::create_quick_operating_point(
            frequency, desiredMagnetizingInductance, ambientTemperature, WaveformLabel::SINUSOIDAL,
            peakToPeak, dutyCycle, dcCurrent);

        json primaryWindingJson = json();
        primaryWindingJson["isolationSide"] = IsolationSide::PRIMARY;
        primaryWindingJson["name"] = "primary";
        primaryWindingJson["numberParallels"] = 1;
        primaryWindingJson["numberTurns"] = numberTurns;
        primaryWindingJson["wire"] = "Dummy";
        OpenMagnetics::Winding primaryCoilFunctionalDescription(primaryWindingJson);
        json CoilFunctionalDescriptionJson = json::array();
        CoilFunctionalDescriptionJson.push_back(primaryWindingJson);
        json windingJson = json();
        windingJson["bobbin"] = "Dummy";
        windingJson["functionalDescription"] = CoilFunctionalDescriptionJson;
        OpenMagnetics::Coil windingAux(windingJson);
        winding = windingAux;

        core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
    }

    TEST_CASE("Test_Inductance_Ferrite_Ground", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 666;
        double frequency = 20000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.003);

        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizingInductanceModel("ZHANG");

        double expectedValue = 23.3e-3;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, -1, gapping, coreShape,
                                coreMaterial, core, winding, inputs);

        auto operatingPoint = inputs.get_operating_point(0);
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(magnetizingInductance, max_error * expectedValue));
    }

    TEST_CASE("Test_Inductance_Ferrite_Web", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        // This tests checks that the operating is not crashing

        json coreData = json::parse(
            R"({"functionalDescription": {"gapping": [{"area": null, "coordinates": null,
            "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 0.001,
            "sectionDimensions": null, "shape": null, "type": "subtractive"}, {"area": null,
            "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null,
            "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}, {"area":
            null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface":
            null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}],
            "material": "3C97", "name": "My Core", "numberStacks": 1, "shape": {"aliases": [],
            "dimensions": {"A": 0.0391, "B": 0.0198, "C": 0.0125, "D": 0.0146, "E": 0.030100000000000002,
            "F": 0.0125, "G": 0.0, "H": 0.0}, "family": "etd", "familySubtype": "1",
            "magneticCircuit": null, "name": "ETD 39/20/13", "type": "standard"}, "type": "two-piece set"}})");
        json windingData =
            json::parse(R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary",
                        "numberParallels": 1, "numberTurns": 1, "wire": "Dummy"}], "layersDescription":
                        null, "sectionsDescription": null, "turnsDescription": null})");
        json operatingPointData = json::parse(
            R"({"conditions": {"ambientRelativeHumidity": null, "ambientTemperature": 25.0, "cooling": null,
            "name": null}, "excitationsPerWinding": [{"current": {"harmonics": null, "processed": null,
            "waveform": {"ancillaryLabel": null, "data": [-5.0, 5.0, -5.0], "numberPeriods": null, "time":
            [0.0, 2.5e-06, 1e-05]}}, "frequency": 100000.0, "magneticField": null, "magneticFluxDensity": null,
            "magnetizingCurrent": null, "name": "My Operating Point", "voltage": {"harmonics": null,
            "processed": null, "waveform": {"ancillaryLabel": null, "data": [7.5, 7.5, -2.5, -2.5, 7.5],
            "numberPeriods": null, "time": [0.0, 2.5e-06, 2.5e-06, 1e-05, 1e-05]}}}],"name": null})");

        Core core(coreData);
        OpenMagnetics::Coil winding(windingData);
        OperatingPoint operatingPoint(operatingPointData);
        MagnetizingInductance magnetizing_inductance("ZHANG");
        magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();
    }

    TEST_CASE("Test_Inductance_Powder_Web", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        // This tests checks that the operating is not crashing

        json coreData = json::parse(
            R"({"functionalDescription": {"gapping": [{"area": null, "coordinates": null,
            "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 0.001,
            "sectionDimensions": null, "shape": null, "type": "subtractive"}, {"area": null,
            "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null,
            "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}, {"area":
            null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface":
            null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}],
            "material": {"bhCycle": null, "curieTemperature": 500.0, "remanence": null, "resistivity": [{"value": 5, "temperature": 20}], "family": "High Flux",
            "manufacturerInfo": {"cost": null, "name": "Magnetics", "reference": null, "status": null},
            "material": "powder", "name": "High Flux 26", "permeability": {"amplitude":
            null, "initial": {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak":
            null, "modifiers": {"EQ/LP": {"frequencyFactor": null, "magneticFieldDcBiasFactor": {"a": 0.01,
            "b": 1.58277e-17, "c": 3.243}, "method": "magnetics", "temperatureFactor": null}, "default":
            {"frequencyFactor": {"a": 0.0, "b": -2.56e-08, "c": 3.4300000000000005e-15, "d": -7.34e-22,
            "e": 3.99e-29}, "magneticFieldDcBiasFactor": {"a": 0.01, "b": 1.02934e-13, "c": 2.426},
            "method": "magnetics", "temperatureFactor": {"a": -0.0033, "b": 0.000129, "c":
            3.799999999999999e-08, "d": 0.0, "e": 0.0}}}, "temperature": null, "tolerance": null, "value":
            26.0}}, "saturation": [{"magneticField": 7957.0, "magneticFluxDensity": 0.9, "temperature":
            100.0}], "type": "commercial", "volumetricLosses": {"EQ/LP": [{"coerciveForce": null,
            "method": "steinmetz", "ranges": [{"alpha": 2.165, "beta": 1.357, "ct0": null, "ct1": null,
            "ct2": null, "k": 14.41908, "maximumFrequency": 1000000000.0, "minimumFrequency": 1.0}],
            "referenceVolumetricLosses": null}], "default":
            [{"coerciveForce": null, "method": "steinmetz", "ranges": [{"alpha": 2.218, "beta": 1.24,
            "ct0": null, "ct1": null, "ct2": null, "k": 93.80774, "maximumFrequency": 1000000000.0,
            "minimumFrequency": 1.0}], "referenceVolumetricLosses": null, "remanence": null, "resistivity":
            null}]}}, "name": "My Core", "numberStacks": 2, "shape": {"aliases": [], "dimensions":
            {"A": 0.0351, "B": 0.0155, "C": 0.01, "D": 0.0095, "E": 0.025, "F": 0.01, "G": 0.0, "H":
            0.0}, "family": "e", "familySubtype": null, "magneticCircuit": "open", "name": "E 35/10",
            "type": "standard"}, "type": "two-piece set"}, "geometricalDescription": null,
            "processedDescription": null})");
        json windingData =
            json::parse(R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary",
                        "numberParallels": 1, "numberTurns": 23, "wire": "Dummy"}], "layersDescription":
                        null, "sectionsDescription": null, "turnsDescription": null})");
        json operatingPointData = json::parse(
            R"({"conditions": {"ambientRelativeHumidity": null, "ambientTemperature": 125.0, "cooling": null,
            "name": null}, "excitationsPerWinding": [{"current": {"harmonics": null, "processed": null,
            "waveform": {"ancillaryLabel": null, "data": [-5.0, 5.0, -5.0], "numberPeriods": null, "time":
            [0.0, 2.5e-06, 1e-05]}}, "frequency": 100000.0, "magneticField": null, "magneticFluxDensity": null,
            "magnetizingCurrent": null, "name": "My Operating Point", "voltage": {"harmonics": null,
            "processed": null, "waveform": {"ancillaryLabel": null, "data": [7.5, 7.5, -2.5, -2.5, 7.5],
            "numberPeriods": null, "time": [0.0, 2.5e-06, 2.5e-06, 1e-05, 1e-05]}}}], "name": null})");

        CoreMaterial coreMaterial(coreData["functionalDescription"]["material"]);
        Core core(coreData);
        OpenMagnetics::Coil winding(windingData);
        OperatingPoint operatingPoint(operatingPointData);
        MagnetizingInductance magnetizing_inductance("ZHANG");
        magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();
    }

    TEST_CASE("Test_Inductance_High_Flux_40_Web", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        // This tests checks that the operating is not crashing

        json coreData = json::parse(R"({"functionalDescription": {"gapping": [{"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 0.001, "sectionDimensions": null, "shape": null, "type": "subtractive"}, {"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}, {"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}], "material": "High Flux 40", "numberStacks": 1, "shape": {"aliases": [], "dimensions": {"A": 0.0391, "B": 0.0198, "C": 0.0125, "D": 0.0146, "E": 0.030100000000000002, "F": 0.0125, "G": 0.0, "H": 0.0}, "family": "etd", "familySubtype": "1", "magneticCircuit": null, "name": "ETD 39/20/13", "type": "standard"}, "type": "two-piece set"}, "geometricalDescription": null, "manufacturerInfo": null, "name": "My Core", "processedDescription": null})");
        json windingData =
            json::parse(R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary", "numberParallels": 1, "numberTurns": 24, "wire": "Dummy"}], "layersDescription": null, "sectionsDescription": null, "turnsDescription": null})");
        json inputsData = json::parse(R"({"designRequirements": {"altitude": null, "cti": null, "insulationType": null, "leakageInductance": null, "magnetizingInductance": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.0001279222825940401}, "name": null, "operatingTemperature": null, "overvoltageCategory": null, "pollutionDegree": null, "turnsRatios": []}, "operatingPoints": [{"conditions": {"ambientRelativeHumidity": null, "ambientTemperature": 25.0, "cooling": null, "name": null}, "excitationsPerWinding": [{"current": {"harmonics": null, "processed": null, "waveform": {"ancillaryLabel": null, "data": [-5.0, 5.0, -5.0], "numberPeriods": null, "time": [0.0, 2.5e-06, 1e-05]}}, "frequency": 100000.0, "magneticFieldStrength": null, "magneticFluxDensity": null, "magnetizingCurrent": null, "name": "My Operating Point", "voltage": null}], "name": null}]})");

        Core core(coreData);
        OpenMagnetics::Coil winding(windingData);
        OpenMagnetics::Inputs inputs(inputsData);
        auto operatingPoint = inputs.get_operating_point(0);
        MagnetizingInductance magnetizing_inductance("ZHANG");
        magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();
    }

    TEST_CASE("Test_Inductance_Ferrite_Spacer", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 666;
        double frequency = 20000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_spacer_gap(0.003);
        double expectedValue = 13.5e-3;

        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, -1, gapping, coreShape,
                                coreMaterial, core, winding, inputs);

        auto operatingPoint = inputs.get_operating_point(0);
        auto aux = magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint);
        double magnetizingInductance = aux.get_magnetizing_inductance().get_nominal().value();

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(magnetizingInductance, max_error * expectedValue));
    }

    TEST_CASE("Test_Inductance_Ferrite_Ground_Few_Turns", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 42;
        double numberTurns = 9;
        double frequency = 100000;
        std::string coreShape = "E 47/20/16";
        std::string coreMaterial = "N87";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0004);
        double expectedValue = 63e-6;

        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, -1, gapping, coreShape,
                                coreMaterial, core, winding, inputs);

        auto operatingPoint = inputs.get_operating_point(0);
        double magnetizingInductance =
            magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(magnetizingInductance, max_error * expectedValue));
    }

    TEST_CASE("Test_Inductance_Powder", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        double dcCurrent = 96;
        double ambientTemperature = 25;
        double numberTurns = 13;
        double frequency = 68000;
        std::string coreShape = "E 42/21/15";
        std::string coreMaterial = "Edge 60";
        auto gapping = OpenMagneticsTesting::get_residual_gap();
        double expectedValue = 15.7e-6;

        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, 20e6, gapping, coreShape,
                                coreMaterial, core, winding, inputs);

        auto operatingPoint = inputs.get_operating_point(0);
        auto aux = magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint);
        double magnetizingInductance = aux.get_magnetizing_inductance().get_nominal().value();

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(magnetizingInductance, max_error * expectedValue));
    }

    TEST_CASE("Test_NumberTurns_Ferrite_Ground", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 23.3e-3;
        double frequency = 20000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.003);

        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedValue = 666;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, -1, desiredMagnetizingInductance, gapping,
                                coreShape, coreMaterial, core, winding, inputs);

        double numberTurns = magnetizing_inductance.calculate_number_turns_from_gapping_and_inductance(core, &inputs);

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(numberTurns, max_error * expectedValue));
    }

    TEST_CASE("Test_NumberTurns_Powder", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        double dcCurrent = 96;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 15.7e-6;
        double frequency = 68000;
        std::string coreShape = "E 42/21/15";
        std::string coreMaterial = "Edge 60";
        auto gapping = OpenMagneticsTesting::get_residual_gap();

        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedValue = 13;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, -1, desiredMagnetizingInductance, gapping,
                                coreShape, coreMaterial, core, winding, inputs);

        double numberTurns = magnetizing_inductance.calculate_number_turns_from_gapping_and_inductance(core, &inputs);

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(numberTurns, max_error * expectedValue));
    }

    TEST_CASE("Test_Gapping_Ferrite_Ground", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 23.3e-3;
        double numberTurns = 666;
        double frequency = 20000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";

        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedValue = 0.003;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, desiredMagnetizingInductance, {},
                                coreShape, coreMaterial, core, winding, inputs);

        auto gapping = magnetizing_inductance.calculate_gapping_from_number_turns_and_inductance(
            core, winding, &inputs, GappingType::GROUND);

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(gapping[0].get_length(), max_error * expectedValue));
    }

    TEST_CASE("Test_Gapping_U_Shape_Ferrite_Ground", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 23.3e-3;
        double numberTurns = 666;
        double frequency = 20000;
        std::string coreShape = "U 26/22/16";
        std::string coreMaterial = "3C97";

        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedValue = 0.0066;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, desiredMagnetizingInductance, {},
                                coreShape, coreMaterial, core, winding, inputs);

        auto gapping = magnetizing_inductance.calculate_gapping_from_number_turns_and_inductance(
            core, winding, &inputs, GappingType::GROUND);

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(gapping[0].get_length(), max_error * expectedValue));
    }

    TEST_CASE("Test_Gapping_Ferrite_Distributed", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 23.3e-3;
        double numberTurns = 666;
        double frequency = 20000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";

        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedValue = 0.0004;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, desiredMagnetizingInductance, {},
                                coreShape, coreMaterial, core, winding, inputs);

        auto gapping = magnetizing_inductance.calculate_gapping_from_number_turns_and_inductance(
            core, winding, &inputs, GappingType::DISTRIBUTED);

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(gapping[0].get_length(), max_error * expectedValue));
        REQUIRE(7UL == gapping.size());
    }

    TEST_CASE("Test_Gapping_Ferrite_Distributed_More_Gap_Precision", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 23.3e-3;
        double numberTurns = 666;
        double frequency = 20000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";

        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedValue = 0.0004;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, desiredMagnetizingInductance, {},
                                coreShape, coreMaterial, core, winding, inputs);

        auto gapping = magnetizing_inductance.calculate_gapping_from_number_turns_and_inductance(
            core, winding, &inputs, GappingType::DISTRIBUTED, 5);

        REQUIRE(expectedValue == gapping[0].get_length());
        REQUIRE(7UL == gapping.size());
    }

    TEST_CASE("Test_Gapping_Classic_Web", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        // This tests checks that the operating is not crashing

        json coreData = json::parse(
            R"({"functionalDescription": {"bobbin": null, "gapping": [{"area": null, "coordinates": null,
            "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 0.001,
            "sectionDimensions": null, "shape": null, "type": "subtractive"}, {"area": null,
            "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null,
            "length": 5e-06, "sectionDimensions": null, "shape": null, "type": "residual"}, {"area":
            null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface":
            null, "length": 5e-06, "sectionDimensions": null, "shape": null, "type": "residual"}],
            "material": "3C95", "name": "My Core", "numberStacks": 1, "shape": {"aliases": ["ETD 54"],
            "dimensions": {"A": 0.0545, "B": 0.0276, "C": 0.0189, "D": 0.0202, "E": 0.0412, "F":
            0.0189, "G": 0.0, "H": 0.0}, "family": "etd", "familySubtype": null, "magneticCircuit":
            "open", "name": "ETD 54/28/19", "type": "standard"}, "type": "two-piece set"},
            "geometricalDescription": [{"coordinates": [0.0, 0.0, 0.0], "dimensions": null,
            "insulationMaterial": null, "machining": [{"coordinates": [0.0, 7.5e-05, 0.0], "length":
            0.00015}, {"coordinates": [0.0, 0.0073, 0.0], "length": 0.0003}], "material": "3C97",
            "rotation": [3.141592653589793, 3.141592653589793, 0.0], "shape": {"aliases": [], "dimensions":
            {"A": 0.0391, "B": 0.0198, "C": 0.0125, "D": 0.0146, "E": 0.030100000000000002, "F": 0.0125},
            "family": "etd", "familySubtype": "1", "magneticCircuit": null, "name": "ETD 39/20/13",
            "type": "standard"}, "type": "half set"}, {"coordinates": [0.0, 0.0, 0.0], "dimensions":
            null, "insulationMaterial": null, "machining": [{"coordinates": [0.0, -0.0073, 0.0], "length":
            0.0003}, {"coordinates": [0.0, -7.5e-05, 0.0], "length": 0.00015}], "material": "3C97",
            "rotation": [0.0, 0.0, 0.0], "shape": {"aliases": [], "dimensions": {"A": 0.0391, "B":
            0.0198, "C": 0.0125, "D": 0.0146, "E": 0.030100000000000002, "F": 0.0125}, "family": "etd",
            "familySubtype": "1", "magneticCircuit": null, "name": "ETD 39/20/13", "type": "standard"},
            "type": "half set"}], "processedDescription": {"columns": [{"area": 0.000123, "coordinates":
            [0.0, 0.0, 0.0], "depth": 0.0125, "height": 0.0292, "shape": "round", "type": "central",
            "width": 0.0125}, {"area": 6.2e-05, "coordinates": [0.017301, 0.0, 0.0], "depth": 0.0125,
            "height": 0.0292, "shape": "irregular", "type": "lateral", "width": 0.004501}, {"area":
            6.2e-05, "coordinates": [-0.017301, 0.0, 0.0], "depth": 0.0125, "height": 0.0292, "shape":
            "irregular", "type": "lateral", "width": 0.004501}], "depth": 0.0125, "effectiveParameters":
            {"effectiveArea": 0.0001249790616277593, "effectiveLength": 0.09385923258669904, "effectiveVolume":
            1.1730438813787252e-05, "minimumArea": 0.0001227184630308513}, "height": 0.0396, "width": 0.0391,
            "windingWindows": [{"angle": null, "area": 0.00025696000000000003, "coordinates": [0.00625, 0.0],
            "height": 0.0292, "radialHeight": null, "width": 0.0088}]}})");
        json windingData =
            json::parse(R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary",
                        "numberParallels": 1, "numberTurns": 40, "wire": "Dummy"}], "layersDescription":
                        null, "sectionsDescription": null, "turnsDescription": null})");
        json inputsData = json::parse(
            R"({"designRequirements": {"altitude": null, "cti": null, "insulationType": null,
            "leakageInductance": null, "magnetizingInductance": {"excludeMaximum": null, "excludeMinimum":
            null, "maximum": null, "minimum": null, "nominal": 0.0004126820555843872}, "name": null,
            "operatingTemperature": null, "overvoltageCategory": null, "pollutionDegree": null,
            "turnsRatios": []}, "operatingPoints": [{"conditions": {"ambientRelativeHumidity": null,
            "ambientTemperature": 25.0, "cooling": null, "name": null}, "excitationsPerWinding":
            [{"current": {"harmonics": null, "processed": null, "waveform": {"ancillaryLabel": null,
            "data": [41.0, 51.0, 41.0], "numberPeriods": null, "time": [0.0, 2.4999999999999998e-06, 1e-05]}},
            "frequency": 100000.0, "magneticField": null, "magneticFluxDensity": null, "magnetizingCurrent":
            null, "name": "My Operating Point", "voltage": {"harmonics": null, "processed": null,
            "waveform": {"ancillaryLabel": null, "data": [7.5, 7.5, -2.4999999999999996, -2.4999999999999996,
            7.5], "numberPeriods": null, "time": [0.0, 2.4999999999999998e-06, 2.4999999999999998e-06, 1e-05,
            1e-05]}}}], "name": null}]})");
        GappingType gappingType =
            magic_enum::enum_cast<GappingType>("DISTRIBUTED").value();

        Core core(coreData);
        OpenMagnetics::Coil winding(windingData);
        OpenMagnetics::Inputs inputs(inputsData);
        MagnetizingInductance magnetizing_inductance("CLASSIC");
        auto gapping =
            magnetizing_inductance.calculate_gapping_from_number_turns_and_inductance(core, winding, &inputs, gappingType, 5);

        REQUIRE(gapping.size() == 5);
    }

    TEST_CASE("Test_Gapping_Web", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        // This tests checks that the operating is not crashing
        json coreData = json::parse(
            R"({"functionalDescription": {"bobbin": null, "gapping": [{"area": 0.000369, "coordinates": [0.0,
            0.00, 0.0], "distanceClosestNormalSurface": 0.022448, "distanceClosestParallelSurface":
            0.011524999999999999, "length": 0.0001, "sectionDimensions": [0.02165, 0.02165], "shape": "round",
            "type": "subtractive"}, {"area": 0.000184, "coordinates": [0.026126, 0.0, 0.0],
            "distanceClosestNormalSurface": 0.022448, "distanceClosestParallelSurface": 0.011524999999999999,
            "length": 5e-06, "sectionDimensions": [0.007551, 0.02165], "shape": "irregular", "type":
            "residual"}, {"area": 0.000184, "coordinates": [-0.026126, 0.0, 0.0],
            "distanceClosestNormalSurface": 0.022448, "distanceClosestParallelSurface": 0.011524999999999999,
            "length": 5e-06, "sectionDimensions": [0.007551, 0.02165], "shape": "irregular", "type":
            "residual"}], "material": "3C95", "name": "My Core", "numberStacks": 1, "shape":
            {"aliases": ["ETD 54"], "dimensions": {"A": 0.0545, "B": 0.0276, "C": 0.0189, "D": 0.0202,
            "E": 0.0412, "F": 0.0189}, "family": "etd", "familySubtype": null, "magneticCircuit":
            "open", "name": "ETD 54/28/19", "type": "standard"}, "type": "two-piece set"},
            "geometricalDescription": null, "processedDescription": null})");
        json windingData =
            json::parse(R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary",
                        "numberParallels": 1, "numberTurns": 1, "wire": "Dummy"}], "layersDescription":
                        null, "sectionsDescription": null, "turnsDescription": null})");
        json inputsData = json::parse(
            R"({"designRequirements": {"altitude": null, "cti": null, "insulationType": null,
            "leakageInductance": null, "magnetizingInductance": {"excludeMaximum": null, "excludeMinimum":
            null, "maximum": null, "minimum": null, "nominal": 0.004654652816558039}, "name": null,
            "operatingTemperature": null, "overvoltageCategory": null, "pollutionDegree": null,
            "turnsRatios": []}, "operatingPoints": [{"conditions": {"ambientRelativeHumidity": null,
            "ambientTemperature": 25.0, "cooling": null, "name": null}, "excitationsPerWinding":
            [{"current": {"harmonics": null, "processed": null, "waveform": {"ancillaryLabel": null,
            "data": [41.0, 51.0, 41.0], "numberPeriods": null, "time": [0.0, 2.4999999999999998e-06, 1e-05]}},
            "frequency": 100000.0, "magneticField": null, "magneticFluxDensity": null, "magnetizingCurrent":
            null, "name": "My Operating Point", "voltage": {"harmonics": null, "processed": null,
            "waveform": {"ancillaryLabel": null, "data": [7.5, 7.5, -2.4999999999999996, -2.4999999999999996,
            7.5], "numberPeriods": null, "time": [0.0, 2.4999999999999998e-06, 2.4999999999999998e-06, 1e-05,
            1e-05]}}}], "name": null}]})");
        GappingType gappingType = magic_enum::enum_cast<GappingType>("GROUND").value();

        Core core(coreData);
        OpenMagnetics::Coil winding(windingData);
        OpenMagnetics::Inputs inputs(inputsData);
        MagnetizingInductance magnetizing_inductance("CLASSIC");
        auto gapping =
            magnetizing_inductance.calculate_gapping_from_number_turns_and_inductance(core, winding, &inputs, gappingType, 5);
    }

    TEST_CASE("Test_Magnetizing_Inductance", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 42;
        double frequency = 20000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);

        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedInductanceValue = 215e-6;
        double currentPeakToPeak = 20;
        double voltagePeakToPeak = 2 * M_PI * frequency * expectedInductanceValue * currentPeakToPeak;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, -1, gapping, coreShape,
                                coreMaterial, core, winding, inputs, voltagePeakToPeak);

        double effectiveArea = core.get_processed_description().value().get_effective_parameters().get_effective_area();
        double expectedMagneticFluxDensity =
            expectedInductanceValue * (currentPeakToPeak / 2) / numberTurns / effectiveArea;

        auto operatingPoint = inputs.get_operating_point(0);
        auto ea = magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint);

        auto magnetizingInductance = ea.first.get_magnetizing_inductance().get_nominal().value();;
        auto magneticFluxDensity = ea.second;

        auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value().get_data();
        auto magneticFluxDensityWaveformPeak =
            *max_element(std::begin(magneticFluxDensityWaveform), std::end(magneticFluxDensityWaveform));
        OperatingPointExcitation primaryExcitation =
            OpenMagnetics::Inputs::get_primary_excitation(operatingPoint);

        REQUIRE_THAT(expectedInductanceValue, Catch::Matchers::WithinAbs(magnetizingInductance, max_error * expectedInductanceValue));
        REQUIRE_THAT(expectedMagneticFluxDensity, Catch::Matchers::WithinAbs(magneticFluxDensityWaveformPeak, max_error * expectedMagneticFluxDensity));
        REQUIRE(bool(primaryExcitation.get_voltage()));
        REQUIRE(bool(primaryExcitation.get_magnetizing_current()));

        if (primaryExcitation.get_current()) {
            auto currentProcessed = primaryExcitation.get_current().value().get_processed().value();
            auto magnetizingCurrentProcessed = primaryExcitation.get_current().value().get_processed().value();
            REQUIRE_THAT(currentPeakToPeak, Catch::Matchers::WithinAbs(operatingPoint.get_mutable_excitations_per_winding()[0]
                            .get_magnetizing_current()
                            .value()
                            .get_processed()
                            .value()
                            .get_peak_to_peak()
                            .value(), max_error * currentPeakToPeak));
        }
    }

    TEST_CASE("Test_Gapping_Web_No_Voltage", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        // This tests checks that the operating is not crashing
        json coreData = json::parse(
            R"({"functionalDescription": {"bobbin": null, "gapping": [{"area": 0.000369, "coordinates": [0.0,
            0.0, 0.0], "distanceClosestNormalSurface": 0.022448, "distanceClosestParallelSurface":
            0.011524999999999999, "length": 0.0001, "sectionDimensions": [0.02165, 0.02165], "shape": "round",
            "type": "subtractive"}, {"area": 0.000184, "coordinates": [0.026126, 0.0, 0.0],
            "distanceClosestNormalSurface": 0.022448, "distanceClosestParallelSurface": 0.011524999999999999,
            "length": 5e-06, "sectionDimensions": [0.007551, 0.02165], "shape": "irregular", "type":
            "residual"}, {"area": 0.000184, "coordinates": [-0.026126, 0.0, 0.0],
            "distanceClosestNormalSurface": 0.022448, "distanceClosestParallelSurface": 0.011524999999999999,
            "length": 5e-06, "sectionDimensions": [0.007551, 0.02165], "shape": "irregular", "type":
            "residual"}], "material": "3C95", "name": "My Core", "numberStacks": 1, "shape":
            {"aliases": ["ETD 54"], "dimensions": {"A": 0.0545, "B": 0.0276, "C": 0.0189, "D": 0.0202,
            "E": 0.0412, "F": 0.0189}, "family": "etd", "familySubtype": null, "magneticCircuit":
            "open", "name": "ETD 54/28/19", "type": "standard"}, "type": "two-piece set"},
            "geometricalDescription": null, "processedDescription": null})");
        json windingData =
            json::parse(R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary",
                        "numberParallels": 1, "numberTurns": 1, "wire": "Dummy"}], "layersDescription":
                        null, "sectionsDescription": null, "turnsDescription": null})");
        json inputsData = json::parse(
            R"({"designRequirements": {"altitude": null, "cti": null, "insulationType": null,
            "leakageInductance": null, "magnetizingInductance": {"excludeMaximum": null, "excludeMinimum":
            null, "maximum": null, "minimum": null, "nominal": 0.00004654652816558039}, "name": null,
            "operatingTemperature": null, "overvoltageCategory": null, "pollutionDegree": null,
            "turnsRatios": []}, "operatingPoints": [{"conditions": {"ambientRelativeHumidity": null,
            "ambientTemperature": 25.0, "cooling": null, "name": null}, "excitationsPerWinding":
            [{"current": {"harmonics": null, "processed": null, "waveform": {"ancillaryLabel": null,
            "data": [41.0, 51.0, 41.0], "numberPeriods": null, "time": [0.0, 2.5e-06, 1e-05]}}, "frequency":
            100000.0, "magneticField": null, "magneticFluxDensity": null, "magnetizingCurrent": null, "name":
            "My Operating Point"}], "name": null}]})");
        GappingType gappingType = magic_enum::enum_cast<GappingType>("GROUND").value();

        Core core(coreData);
        OpenMagnetics::Coil winding(windingData);
        OpenMagnetics::Inputs inputs(inputsData);
        MagnetizingInductance magnetizing_inductance("CLASSIC");
        auto gapping =
            magnetizing_inductance.calculate_gapping_from_number_turns_and_inductance(core, winding, &inputs, gappingType, 5);
        auto primaryExcitation = inputs.get_operating_point(0).get_mutable_excitations_per_winding()[0];
        double currentPeakToPeak = 10;

        REQUIRE(bool(primaryExcitation.get_voltage()));
        REQUIRE(bool(primaryExcitation.get_current()));
        REQUIRE(bool(primaryExcitation.get_magnetizing_current()));

        if (primaryExcitation.get_current()) {
            auto currentProcessed = primaryExcitation.get_current().value().get_processed().value();
            auto magnetizingCurrentProcessed = primaryExcitation.get_current().value().get_processed().value();
            REQUIRE_THAT(currentPeakToPeak, Catch::Matchers::WithinAbs(inputs.get_operating_point(0)
                            .get_mutable_excitations_per_winding()[0]
                            .get_magnetizing_current()
                            .value()
                            .get_processed()
                            .value()
                            .get_peak_to_peak()
                            .value(), max_error * currentPeakToPeak));
            REQUIRE_THAT(currentPeakToPeak, Catch::Matchers::WithinAbs(inputs.get_operating_point(0)
                            .get_mutable_excitations_per_winding()[0]
                            .get_current()
                            .value()
                            .get_processed()
                            .value()
                            .get_peak_to_peak()
                            .value(), max_error * currentPeakToPeak));
        }
    }

    TEST_CASE("Test_Inductance_Ferrite_Web_No_Voltage", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        // This tests checks that the operating is not crashing

        json coreData = json::parse(
            R"({"functionalDescription": {"gapping": [{"area": null, "coordinates": null,
            "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 0.001,
            "sectionDimensions": null, "shape": null, "type": "subtractive"}, {"area": null,
            "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null,
            "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}, {"area":
            null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface":
            null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}],
            "material": "3C97", "name": "My Core", "numberStacks": 1, "shape": {"aliases": [],
            "dimensions": {"A": 0.0391, "B": 0.0198, "C": 0.0125, "D": 0.0146, "E": 0.030100000000000002,
            "F": 0.0125, "G": 0.0, "H": 0.0}, "family": "etd", "familySubtype": "1",
            "magneticCircuit": null, "name": "ETD 39/20/13", "type": "standard"}, "type": "two-piece set"}})");
        json windingData =
            json::parse(R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary",
                        "numberParallels": 1, "numberTurns": 10, "wire": "Dummy"}], "layersDescription":
                        null, "sectionsDescription": null, "turnsDescription": null})");
        json operatingPointData = json::parse(
            R"({"conditions": {"ambientRelativeHumidity": null, "ambientTemperature": 25.0, "cooling": null,
            "name": null}, "excitationsPerWinding": [{"current": {"harmonics": null, "processed": null,
            "waveform": {"ancillaryLabel": null, "data": [-5.0, 5.0, -5.0], "numberPeriods": null, "time":
            [0.0, 2.5e-06, 1e-05]}}, "frequency": 100000.0, "magneticField": null, "magneticFluxDensity": null,
            "magnetizingCurrent": null, "name": "My Operating Point"}],"name": null})");

        Core core(coreData);
        OpenMagnetics::Coil winding(windingData);
        OperatingPoint operatingPoint(operatingPointData);
        MagnetizingInductance magnetizing_inductance("ZHANG");
        magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();
        auto primaryExcitation = operatingPoint.get_mutable_excitations_per_winding()[0];
        double currentPeakToPeak = 10;
        double voltagePeakToPeak = 105;

        REQUIRE(bool(primaryExcitation.get_voltage()));
        REQUIRE(bool(primaryExcitation.get_current()));
        REQUIRE(bool(primaryExcitation.get_magnetizing_current()));

        if (primaryExcitation.get_current()) {
            auto currentProcessed = primaryExcitation.get_current().value().get_processed().value();
            auto magnetizingCurrentProcessed = primaryExcitation.get_current().value().get_processed().value();
            REQUIRE_THAT(currentPeakToPeak, Catch::Matchers::WithinAbs(operatingPoint.get_mutable_excitations_per_winding()[0]
                            .get_magnetizing_current()
                            .value()
                            .get_processed()
                            .value()
                            .get_peak_to_peak()
                            .value(), max_error * currentPeakToPeak));
            REQUIRE_THAT(currentPeakToPeak, Catch::Matchers::WithinAbs(operatingPoint.get_mutable_excitations_per_winding()[0]
                            .get_current()
                            .value()
                            .get_processed()
                            .value()
                            .get_peak_to_peak()
                            .value(), max_error * currentPeakToPeak));
        }
        if (primaryExcitation.get_voltage()) {
            REQUIRE_THAT(voltagePeakToPeak, Catch::Matchers::WithinAbs(operatingPoint.get_mutable_excitations_per_winding()[0]
                            .get_voltage()
                            .value()
                            .get_processed()
                            .value()
                            .get_peak_to_peak()
                            .value(), max_error * voltagePeakToPeak));
        }
    }

    TEST_CASE("Test_Magnetizing_Inductance_Toroid", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        settings.reset();
        clear_databases();
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 42;
        double frequency = 20000;
        std::string coreShape = "T 58/41/18";
        std::string coreMaterial = "3C95";
        std::vector<CoreGap> gapping = {};

        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedValue = 6.6e-3;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, -1, gapping, coreShape,
                                coreMaterial, core, winding, inputs);

        auto operatingPoint = inputs.get_operating_point(0);
        double magnetizingInductance =
            magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(magnetizingInductance, max_error * expectedValue));
    }

    TEST_CASE("Test_Magnetizing_Inductance_Toroid_Stacks", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        settings.reset();
        clear_databases();
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 42;
        double frequency = 20000;
        std::string coreShape = "T 58/41/18";
        std::string coreMaterial = "3C95";
        std::vector<CoreGap> gapping = {};

        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedValue = 6.6e-3;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, -1, gapping, coreShape,
                                coreMaterial, core, winding, inputs);

        auto operatingPoint = inputs.get_operating_point(0);
        double magnetizingInductance =
            magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(magnetizingInductance, max_error * expectedValue));

        core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, 2, coreMaterial);
        double magnetizingInductance2Stacks =
            magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();

        expectedValue = magnetizingInductance * 2;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(magnetizingInductance2Stacks, max_error * expectedValue));
    }

    TEST_CASE("Test_Magnetizing_Inductance_RM14_20", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 29;
        double frequency = 100000;
        std::string coreShape = "RM 14/20";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);

        Core core = json::parse(R"({"name": "My Core", "functionalDescription": {"type": "two-piece set", "material": "3C97", "shape": {"aliases": ["RM 14LP", "RM 14/ILP", "RM 14/LP"], "dimensions": {"A": {"minimum": 0.0408, "maximum": 0.0422 }, "B": {"minimum": 0.010150000000000001, "maximum": 0.01025 }, "C": {"minimum": 0.018400000000000003, "maximum": 0.019000000000000003 }, "D": {"minimum": 0.00555, "maximum": 0.00585 }, "E": {"minimum": 0.029, "maximum": 0.0302 }, "F": {"minimum": 0.014400000000000001, "maximum": 0.015000000000000001 }, "G": {"minimum": 0.017 }, "H": {"minimum": 0.0054, "maximum": 0.005600000000000001 }, "J": {"minimum": 0.0335, "maximum": 0.0347 }, "R": {"maximum": 0.00030000000000000003 } }, "family": "rm", "familySubtype": "3", "name": "RM 14/20", "type": "standard", "magneticCircuit": "open"}, "gapping": [{"type": "subtractive", "length": 0.001 }, {"length": 0.000005, "type": "residual"}, {"length": 0.000005, "type": "residual"}, {"length": 0.000005, "type": "residual"} ], "numberStacks": 1 }, "geometricalDescription": null, "processedDescription": null })");
        OpenMagnetics::Coil winding; 
        OpenMagnetics::Inputs inputs; 
        MagnetizingInductance magnetizing_inductance("ZHANG");

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, -1, gapping, coreShape,
                                coreMaterial, core, winding, inputs);

        auto operatingPoint = inputs.get_operating_point(0);
        magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();
    }

    TEST_CASE("Test_Magnetizing_Inductance_Error_Web_0", "[physical-model][magnetizing-inductance][bug][smoke-test]") {
        settings.reset();
        clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 10;
        double frequency = 20000;
        std::string coreShape = "E 65/32/27";
        std::string coreMaterial = "N95";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);

        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedValue = 19e-6;

        int numberStacks = 2;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, -1, gapping, coreShape,
                                coreMaterial, core, winding, inputs, 20, numberStacks);

        auto operatingPoint = inputs.get_operating_point(0);
        double magnetizingInductance =
            magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(magnetizingInductance, max_error * expectedValue));
    }

    TEST_CASE("Test_Magnetizing_Inductance_Error_Web_1", "[physical-model][magnetizing-inductance][bug][smoke-test]") {
        settings.reset();
        clear_databases();

        Core core = json::parse(R"({"name": "650-4637", "functionalDescription": {"type": "two-piece set", "material": "TP4A", "shape": {"aliases": ["E 16/5", "EF 16"], "dimensions": {"A": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0167, "minimum": 0.0155, "nominal": null}, "B": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0082, "minimum": 0.0079, "nominal": null}, "C": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0047, "minimum": 0.0043, "nominal": null}, "D": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0061, "minimum": 0.0057, "nominal": null}, "E": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0119, "minimum": 0.0113, "nominal": null}, "F": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0047, "minimum": 0.0044, "nominal": null}}, "family": "e", "familySubtype": null, "magneticCircuit": "open", "name": "E 16/8/5", "type": "standard"}, "gapping": [], "numberStacks": 1}, "processedDescription": {"columns": [{"area": 2.1e-05, "coordinates": [0.0, 0.0, 0.0], "depth": 0.004501, "height": 0.011802, "minimumDepth": null, "minimumWidth": null, "shape": "rectangular", "type": "central", "width": 0.00455}, {"area": 1.1e-05, "coordinates": [0.006925, 0.0, 0.0], "depth": 0.004501, "height": 0.011802, "minimumDepth": null, "minimumWidth": null, "shape": "rectangular", "type": "lateral", "width": 0.002251}, {"area": 1.1e-05, "coordinates": [-0.006925, 0.0, 0.0], "depth": 0.004501, "height": 0.011802, "minimumDepth": null, "minimumWidth": null, "shape": "rectangular", "type": "lateral", "width": 0.002251}], "depth": 0.0045000000000000005, "effectiveParameters": {"effectiveArea": 2.0062091987236854e-05, "effectiveLength": 0.03756497447228765, "effectiveVolume": 7.53631973361239e-07, "minimumArea": 1.935000000000001e-05}, "height": 0.016100000000000003, "width": 0.0161, "windingWindows": [{"angle": null, "area": 4.1595e-05, "coordinates": [0.002275, 0.0], "height": 0.011800000000000001, "radialHeight": null, "sectionsAlignment": null, "sectionsOrientation": null, "shape": null, "width": 0.0035249999999999995}]}})");
        OpenMagnetics::Coil coil = json::parse(R"({"bobbin": "Dummy", "functionalDescription": [{"name": "PRI", "numberTurns": 192, "numberParallels": 1, "connections": [{"pinName": "2"}, {"pinName": "1"}], "isolationSide": "primary", "wire": "Round 35.0 - Heavy Build"}, {"name": "SEC", "numberTurns": 36, "numberParallels": 1, "connections": [{"pinName": "8"}, {"pinName": "7"}], "isolationSide": "secondary", "wire": "Round 29.0 - Single Build"}, {"name": "AUX", "numberTurns": 20, "numberParallels": 1, "connections": [{"pinName": "4"}, {"pinName": "3"}], "isolationSide": "tertiary", "wire": "Round 35.0 - Heavy Build"}]})");
        OpenMagnetics::Inputs inputs = json::parse(R"({"designRequirements": {"name": "basicRequirements", "magnetizingInductance": {"nominal": 0.00232}, "turnsRatios": [{"nominal": 0.1875}, {"nominal": 0.10416666666666667}]}, "operatingPoints": []})");
        json modelsData = json::parse("{}");

        std::map<std::string, std::string> models = modelsData.get<std::map<std::string, std::string>>();
        GappingType gappingType = magic_enum::enum_cast<GappingType>("GROUND").value();
        
        auto reluctanceModelName = Defaults().reluctanceModelDefault;
        if (models.find("reluctance") != models.end()) {
            std::string modelNameJsonUpper = models["reluctance"];
            std::transform(modelNameJsonUpper.begin(), modelNameJsonUpper.end(), modelNameJsonUpper.begin(), ::toupper);
            reluctanceModelName = magic_enum::enum_cast<ReluctanceModels>(modelNameJsonUpper).value();
        }

        MagnetizingInductance magnetizingInductanceObj(reluctanceModelName);
        std::vector<CoreGap> gapping = magnetizingInductanceObj.calculate_gapping_from_number_turns_and_inductance(core,
                                                                                                           coil,
                                                                                                           &inputs,
                                                                                                           gappingType,
                                                                                                           6);

        core.set_processed_description(std::nullopt);
        core.set_geometrical_description(std::nullopt);
        core.get_mutable_functional_description().set_gapping(gapping);
        core.process_data();
        core.process_gap();
        auto geometricalDescription = core.create_geometrical_description();
        core.set_geometrical_description(geometricalDescription);

        json result;
        to_json(result, core);
    }

    TEST_CASE("Test_Inductance_Powder_E_65", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        double max_error = 0.15;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 10;
        double frequency = 100000;
        std::string coreShape = "E 65/32/27";
        std::string coreMaterial = "Kool Mµ 40";
        auto gapping = OpenMagneticsTesting::get_residual_gap();
        double expectedValue = 23e-6;

        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, 20e6, gapping, coreShape,
                                coreMaterial, core, winding, inputs);

        settings.set_magnetizing_inductance_include_air_inductance(true);

        auto operatingPoint = inputs.get_operating_point(0);
        auto aux = magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint);
        double magnetizingInductance = aux.get_magnetizing_inductance().get_nominal().value();

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(magnetizingInductance, max_error * expectedValue));
        settings.reset();
    }

    TEST_CASE("Test_Inductance_Powder_E_34", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 10;
        double frequency = 100000;
        std::string coreShape = "E 34/14/9";
        std::string coreMaterial = "Edge 26";
        auto gapping = OpenMagneticsTesting::get_residual_gap();
        double expectedValue = 5.6e-6;

        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, 20e6, gapping, coreShape,
                                coreMaterial, core, winding, inputs);

        settings.set_magnetizing_inductance_include_air_inductance(true);

        auto operatingPoint = inputs.get_operating_point(0);
        auto aux = magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint);
        double magnetizingInductance = aux.get_magnetizing_inductance().get_nominal().value();

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(magnetizingInductance, max_error * expectedValue));
        settings.reset();
    }

    TEST_CASE("Test_Inductance_Bug_Web_0", "[physical-model][magnetizing-inductance][bug][smoke-test]") {
        settings.reset();
        clear_databases();

        // This tests checks that the operating is not crashing

        json coreData = json::parse(R"({"distributorsInfo":null,"functionalDescription":{"coating":null,"gapping":[],"material":"MPP 26","numberStacks":1,"shape":{"aliases":["R 80/20/50"],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.08},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.05},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.02}},"family":"t","familySubtype":null,"magneticCircuit":"closed","name":"T 80/20/50","type":"standard"},"type":"toroidal"},"geometricalDescription":[{"coordinates":[0.0,0.0,0.0],"dimensions":null,"insulationMaterial":null,"machining":null,"material":"MPP 26","rotation":[1.5707963267948966,1.5707963267948966,0.0],"shape":{"aliases":["R 80/20/50"],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.08},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.05},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.02}},"family":"t","familySubtype":null,"magneticCircuit":"closed","name":"T 80/20/50","type":"standard"},"type":"toroidal"}],"manufacturerInfo":null,"name":"Custom","processedDescription":{"columns":[{"area":0.0003,"coordinates":[0.0,0.0,0.0],"depth":0.02,"height":0.20420352248333656,"minimumDepth":null,"minimumWidth":null,"shape":"rectangular","type":"central","width":0.015}],"depth":0.02,"effectiveParameters":{"effectiveArea":0.0003,"effectiveLength":0.20420352248333654,"effectiveVolume":6.126105674500096e-05,"minimumArea":0.0003},"height":0.08,"thermalResistance":null,"width":0.08,"windingWindows":[{"angle":360.0,"area":0.001963495408493621,"coordinates":[0.015,0.0],"height":null,"radialHeight":0.025,"sectionsAlignment":null,"sectionsOrientation":null,"shape":null,"width":null}]}})");
        json windingData = json::parse(R"({"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.01,"columnShape":"rectangular","columnThickness":0.0,"columnWidth":0.0075,"coordinates":[0.0,0.0,0.0],"pins":null,"wallThickness":0.0,"windingWindows":[{"angle":360.0,"area":0.001963495408493621,"coordinates":[0.025,0.0,0.0],"height":null,"radialHeight":0.025,"sectionsAlignment":"inner or top","sectionsOrientation":"overlapping","shape":"round","width":null}]}},"functionalDescription":[{"connections":null,"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":1,"wire":{"coating":{"breakdownVoltage":4600.0,"grade":2,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1.8834326265752323e-05},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.004897},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Nearson","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 4.5 - Heavy Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.004989},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"4.5 AWG","strand":null,"type":"round"}},{"connections":null,"isolationSide":"primary","name":"Secondary","numberParallels":1,"numberTurns":1,"wire":{"coating":{"breakdownVoltage":4600.0,"grade":2,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1.8834326265752323e-05},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.004897},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Nearson","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 4.5 - Heavy Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.004989},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"4.5 AWG","strand":null,"type":"round"}},{"connections":null,"isolationSide":"primary","name":"Tertiary","numberParallels":1,"numberTurns":1,"wire":{"coating":{"breakdownVoltage":4600.0,"grade":2,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1.8834326265752323e-05},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.004897},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Nearson","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 4.5 - Heavy Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.004989},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"4.5 AWG","strand":null,"type":"round"}}],"layersDescription":[{"additionalCoordinates":null,"coordinateSystem":"polar","coordinates":[0.0024945,180.0],"dimensions":[0.004989,360.0],"fillingFactor":0.03967937689697995,"insulationMaterial":null,"name":"Primary section 0 layer 0","orientation":"overlapping","partialWindings":[{"connections":null,"parallelsProportion":[1.0],"winding":"Primary"}],"section":"Primary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveParallels"},{"additionalCoordinates":[[-0.0200015,180.0]],"coordinateSystem":"polar","coordinates":[0.0050015,180.0],"dimensions":[2.5e-05,360.0],"fillingFactor":1.0,"insulationMaterial":null,"name":"Insulation between Primary and Primary section 1 layer 0","orientation":"overlapping","partialWindings":[],"section":"Insulation between Primary and Primary section 1","turnsAlignment":"spread","type":"insulation","windingStyle":null},{"additionalCoordinates":null,"coordinateSystem":"polar","coordinates":[0.007508500000000001,180.0],"dimensions":[0.004989,360.0],"fillingFactor":0.03967937689697995,"insulationMaterial":null,"name":"Secondary section 0 layer 0","orientation":"overlapping","partialWindings":[{"connections":null,"parallelsProportion":[1.0],"winding":"Secondary"}],"section":"Secondary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveParallels"},{"additionalCoordinates":[[-0.0250155,180.0]],"coordinateSystem":"polar","coordinates":[0.0100155,180.0],"dimensions":[2.5e-05,360.0],"fillingFactor":1.0,"insulationMaterial":null,"name":"Insulation between Secondary and Secondary section 3 layer 0","orientation":"overlapping","partialWindings":[],"section":"Insulation between Secondary and Secondary section 3","turnsAlignment":"spread","type":"insulation","windingStyle":null},{"additionalCoordinates":null,"coordinateSystem":"polar","coordinates":[0.0125225,180.0],"dimensions":[0.004989,360.0],"fillingFactor":0.03967937689697995,"insulationMaterial":null,"name":"Tertiary section 0 layer 0","orientation":"overlapping","partialWindings":[{"connections":null,"parallelsProportion":[1.0],"winding":"Tertiary"}],"section":"Tertiary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveParallels"},{"additionalCoordinates":[[-0.0300295,180.0]],"coordinateSystem":"polar","coordinates":[0.0150295,180.0],"dimensions":[2.5e-05,360.0],"fillingFactor":1.0,"insulationMaterial":null,"name":"Insulation between Tertiary and Tertiary section 5 layer 0","orientation":"overlapping","partialWindings":[],"section":"Insulation between Tertiary and Tertiary section 5","turnsAlignment":"spread","type":"insulation","windingStyle":null}],"sectionsDescription":[{"coordinateSystem":"polar","coordinates":[0.0024945,180.0],"dimensions":[0.004989,360.0],"fillingFactor":0.035281331722710724,"layersAlignment":null,"layersOrientation":"overlapping","margin":[0.0,0.0],"name":"Primary section 0","partialWindings":[{"connections":null,"parallelsProportion":[1.0],"winding":"Primary"}],"type":"conduction","windingStyle":"windByConsecutiveParallels"},{"coordinateSystem":"polar","coordinates":[0.005001500000000001,180.0],"dimensions":[2.5e-05,360.0],"fillingFactor":1.0,"layersAlignment":null,"layersOrientation":"overlapping","margin":null,"name":"Insulation between Primary and Primary section 1","partialWindings":[],"type":"insulation","windingStyle":null},{"coordinateSystem":"polar","coordinates":[0.007508500000000001,180.0],"dimensions":[0.004989,360.0],"fillingFactor":0.04539484956038453,"layersAlignment":null,"layersOrientation":"overlapping","margin":[0.0,0.0],"name":"Secondary section 0","partialWindings":[{"connections":null,"parallelsProportion":[1.0],"winding":"Secondary"}],"type":"conduction","windingStyle":"windByConsecutiveParallels"},{"coordinateSystem":"polar","coordinates":[0.010015500000000002,180.0],"dimensions":[2.5e-05,360.0],"fillingFactor":1.0,"layersAlignment":null,"layersOrientation":"overlapping","margin":null,"name":"Insulation between Secondary and Secondary section 3","partialWindings":[],"type":"insulation","windingStyle":null},{"coordinateSystem":"polar","coordinates":[0.0125225,180.0],"dimensions":[0.004989,360.0],"fillingFactor":0.06363646652658514,"layersAlignment":null,"layersOrientation":"overlapping","margin":[0.0,0.0],"name":"Tertiary section 0","partialWindings":[{"connections":null,"parallelsProportion":[1.0],"winding":"Tertiary"}],"type":"conduction","windingStyle":"windByConsecutiveParallels"},{"coordinateSystem":"polar","coordinates":[0.015029500000000003,180.0],"dimensions":[2.5e-05,360.0],"fillingFactor":1.0,"layersAlignment":null,"layersOrientation":"overlapping","margin":null,"name":"Insulation between Tertiary and Tertiary section 5","partialWindings":[],"type":"insulation","windingStyle":null}],"turnsDescription":[{"additionalCoordinates":[[-0.042494500000000004,5.204075340636721e-18]],"angle":null,"coordinateSystem":"cartesian","coordinates":[-0.0225055,2.756128853821076e-18],"dimensions":[0.004989,0.004989],"layer":"Primary section 0 layer 0","length":0.08567340574875948,"name":"Primary parallel 0 turn 0","orientation":"clockwise","parallel":0,"rotation":180.0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":[[-0.0475085,5.818113245729203e-18]],"angle":null,"coordinateSystem":"cartesian","coordinates":[-0.0174915,2.142090948728593e-18],"dimensions":[0.004989,0.004989],"layer":"Secondary section 0 layer 0","length":0.11717729687895795,"name":"Secondary parallel 0 turn 0","orientation":"clockwise","parallel":0,"rotation":180.0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":[[-0.0525225,6.432151150821686e-18]],"angle":null,"coordinateSystem":"cartesian","coordinates":[-0.0124775,1.52805304363611e-18],"dimensions":[0.004989,0.004989],"layer":"Tertiary section 0 layer 0","length":0.14868118800915636,"name":"Tertiary parallel 0 turn 0","orientation":"clockwise","parallel":0,"rotation":180.0,"section":"Tertiary section 0","winding":"Tertiary"}]} )");
        std::string file_path_1009 = std::source_location::current().file_name();
        auto json_path_1009 = file_path_1009.substr(0, file_path_1009.rfind("/")).append("/testData/test_inductance_bug_web_0_1009.json");
        std::ifstream json_file_1009(json_path_1009);
        json operatingPointData = json::parse(json_file_1009);

        Core core(coreData);
        OpenMagnetics::Coil winding(windingData);
        OperatingPoint operatingPoint(operatingPointData);
        MagnetizingInductance magnetizingInductance("ZHANG");
        magnetizingInductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();
    }

    TEST_CASE("Test_Inductance_Bug_Web_1", "[physical-model][magnetizing-inductance][bug]") {
        settings.reset();
        clear_databases();

        json coreData = json::parse(R"({"distributorsInfo":null,"functionalDescription":{"coating":null,"gapping":[{"area":0.000315,"coordinates":[0,0,0],"distanceClosestNormalSurface":0.014498,"distanceClosestParallelSurface":0.011999999999999999,"length":0.000005,"sectionDimensions":[0.02,0.02],"shape":"round","type":"residual"},{"area":0.000164,"coordinates":[0.024563,0,0],"distanceClosestNormalSurface":0.014498,"distanceClosestParallelSurface":0.011999999999999999,"length":0.000005,"sectionDimensions":[0.005125,0.032],"shape":"irregular","type":"residual"},{"area":0.000164,"coordinates":[-0.024563,0,0],"distanceClosestNormalSurface":0.014498,"distanceClosestParallelSurface":0.011999999999999999,"length":0.000005,"sectionDimensions":[0.005125,0.032],"shape":"irregular","type":"residual"}],"material":"3C95","numberStacks":1,"shape":{"aliases":["EQ 50/32/20.0","EQ 50/20/32"],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.05},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.02},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.032},"D":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0145},"E":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.044},"F":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.02},"G":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.03205}},"family":"eq","familySubtype":null,"magneticCircuit":"open","name":"EQ 50/32/20","type":"standard"},"type":"two-piece set","magneticCircuit":"open"},"geometricalDescription":[{"coordinates":[0,0,0],"dimensions":null,"insulationMaterial":null,"machining":null,"material":"3C95","rotation":[3.141592653589793,3.141592653589793,0],"shape":{"aliases":["EQ 50/32/20.0","EQ 50/20/32"],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.05},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.02},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.032},"D":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0145},"E":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.044},"F":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.02},"G":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.03205}},"family":"eq","familySubtype":null,"magneticCircuit":"open","name":"EQ 50/32/20","type":"standard"},"type":"half set"},{"coordinates":[0,0,0],"dimensions":null,"insulationMaterial":null,"machining":null,"material":"3C95","rotation":[0,0,0],"shape":{"aliases":["EQ 50/32/20.0","EQ 50/20/32"],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.05},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.02},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.032},"D":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0145},"E":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.044},"F":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.02},"G":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.03205}},"family":"eq","familySubtype":null,"magneticCircuit":"open","name":"EQ 50/32/20","type":"standard"},"type":"half set"}],"manufacturerInfo":null,"name":"Custom","processedDescription":{"columns":[{"area":0.000315,"coordinates":[0,0,0],"depth":0.02,"height":0.029,"minimumDepth":null,"minimumWidth":null,"shape":"round","type":"central","width":0.02},{"area":0.000164,"coordinates":[0.024563,0,0],"depth":0.032,"height":0.029,"minimumDepth":null,"minimumWidth":0.003001,"shape":"irregular","type":"lateral","width":0.005125},{"area":0.000164,"coordinates":[-0.024563,0,0],"depth":0.032,"height":0.029,"minimumDepth":null,"minimumWidth":0.003001,"shape":"irregular","type":"lateral","width":0.005125}],"depth":0.032,"effectiveParameters":{"effectiveArea":0.0003298035730425377,"effectiveLength":0.10383305467139949,"effectiveVolume":0.00003424451243054872,"minimumArea":0.0003141592653589793},"height":0.04,"thermalResistance":null,"width":0.05,"windingWindows":[{"angle":null,"area":0.00034799999999999995,"coordinates":[0.01,0],"height":0.029,"radialHeight":null,"sectionsAlignment":null,"sectionsOrientation":null,"shape":null,"width":0.011999999999999999}]}})");
        json windingData = json::parse(R"({"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.012082738036576609,"columnShape":"round","columnThickness":0.0020827380365766087,"columnWidth":0.012082738036576609,"coordinates":[0,0,0],"pins":null,"wallThickness":0.0017356150304805081,"windingWindows":[{"angle":null,"area":0.00025317549908941776,"coordinates":[0.017041369018288302,0,0],"height":0.025528769939038985,"radialHeight":null,"sectionsAlignment":null,"sectionsOrientation":null,"shape":"rectangular","width":0.00991726196342339}]}},"functionalDescription":[{"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":1,"wire":"Dummy"}],"layersDescription":null,"sectionsDescription":null,"turnsDescription":null})");
        json operatingPointData = json::parse(R"({"name":"Op. Point No. 1","conditions":{"ambientTemperature":100},"excitationsPerWinding":[{"frequency":50,"current":{"harmonics":{"amplitudes":[0,31.21831,0.348898,0.681642,11.946624,12.028013,10.143728,9.443625,10.985242,10.23761,10.161232,10.519392,11.473431,12.334118,15.790292,12.339356,3.782311,3.522273,3.435625,3.592056,3.338095,0.365192],"frequencies":[0,50,15000,25000,35000,45000,55000,65000,75000,85000,95000,105000,115000,125000,135000,145000,155000,165000,175000,185000,195000,205000]}},"voltage":{"harmonics":{"amplitudes":[0,1.136963,1.415776,2.520523,36.863994,43.911459,45.631752,50.314416,67.374707,71.317378,79.215135,90.682108,108.339425,126.882128,176.340424,143.759392,47.845538,47.502176,49.241268,54.317715,53.222847,6.05556],"frequencies":[0,50,15000,25000,35000,45000,55000,65000,75000,85000,95000,105000,115000,125000,135000,145000,155000,165000,175000,185000,195000,205000]}}}]})");

        Core core(coreData);
        OpenMagnetics::Coil winding(windingData);
        OperatingPoint operatingPoint(operatingPointData);
        MagnetizingInductance magnetizingInductance("ZHANG");
        magnetizingInductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();
    }

    TEST_CASE("Test_Magnetizing_Inductance_Error_Web_2", "[physical-model][magnetizing-inductance][bug][smoke-test]") {
        settings.reset();
        clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 16;
        double frequency = 570000;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "DMR51W";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.0005, 4);

        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedValue = 30e-6;

        int numberStacks = 1;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, -1, gapping, coreShape,
                                coreMaterial, core, winding, inputs, 20, numberStacks);

        auto operatingPoint = inputs.get_operating_point(0);
        double magnetizingInductance =
            magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(magnetizingInductance, max_error * expectedValue));
    }

    TEST_CASE("Test_Magnetizing_Inductance_Error_Web_3", "[physical-model][magnetizing-inductance][bug][smoke-test]") {
        settings.reset();
        clear_databases();

        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "DMR51W";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.0005, 4);

        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");
        double expectedValue = 30e-6;
        std::string file_path = std::source_location::current().file_name();

        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/Error_inductance_with_Csv.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        auto operatingPoint = mas.get_inputs().get_operating_points()[0];

        double magnetizingInductance =
            magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(mas.get_magnetic().get_core(), mas.get_magnetic().get_coil(), &operatingPoint).get_magnetizing_inductance().get_nominal().value();

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(magnetizingInductance, max_error * expectedValue));
    }

}  // namespace
