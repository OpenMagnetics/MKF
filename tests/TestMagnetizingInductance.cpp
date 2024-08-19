#include "MagnetizingInductance.h"
#include "BobbinWrapper.h"
#include "TestingUtils.h"
#include "Settings.h"
#include "Utils.h"
#include "json.hpp"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <typeinfo>
#include <vector>

using json = nlohmann::json;

SUITE(MagnetizingInductance) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    double max_error = 0.05;
    void prepare_test_parameters(double dcCurrent, double ambientTemperature, double frequency, double numberTurns,
                                 double desiredMagnetizingInductance, std::vector<OpenMagnetics::CoreGap> gapping,
                                 std::string coreShape, std::string coreMaterial, OpenMagnetics::CoreWrapper& core,
                                 OpenMagnetics::CoilWrapper& winding, OpenMagnetics::InputsWrapper& inputs,
                                 double peakToPeak = 20, int numberStacks = 1) {
        double dutyCycle = 0.5;

        inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(
            frequency, desiredMagnetizingInductance, ambientTemperature, OpenMagnetics::WaveformLabel::SINUSOIDAL,
            peakToPeak, dutyCycle, dcCurrent);

        json primaryWindingJson = json();
        primaryWindingJson["isolationSide"] = OpenMagnetics::IsolationSide::PRIMARY;
        primaryWindingJson["name"] = "primary";
        primaryWindingJson["numberParallels"] = 1;
        primaryWindingJson["numberTurns"] = numberTurns;
        primaryWindingJson["wire"] = "Dummy";
        OpenMagnetics::CoilFunctionalDescription primaryCoilFunctionalDescription(primaryWindingJson);
        json CoilFunctionalDescriptionJson = json::array();
        CoilFunctionalDescriptionJson.push_back(primaryWindingJson);
        json windingJson = json();
        windingJson["bobbin"] = "Dummy";
        windingJson["functionalDescription"] = CoilFunctionalDescriptionJson;
        OpenMagnetics::CoilWrapper windingAux(windingJson);
        winding = windingAux;

        core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
    }

    TEST(Test_Inductance_Ferrite_Ground) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 666;
        double frequency = 20000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.003);

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::CoilWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedValue = 23.3e-3;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, -1, gapping, coreShape,
                                coreMaterial, core, winding, inputs);

        auto operatingPoint = inputs.get_operating_point(0);
        double magnetizingInductance =
            magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();

        CHECK_CLOSE(expectedValue, magnetizingInductance, max_error * expectedValue);
    }

    TEST(Test_Inductance_Ferrite_Web) {
        settings->reset();
        OpenMagnetics::clear_databases();

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

        OpenMagnetics::CoreWrapper core(coreData);
        OpenMagnetics::CoilWrapper winding(windingData);
        OpenMagnetics::OperatingPoint operatingPoint(operatingPointData);
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");
        magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();
    }

    TEST(Test_Inductance_Powder_Web) {
        settings->reset();
        OpenMagnetics::clear_databases();

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

        OpenMagnetics::CoreMaterial coreMaterial(coreData["functionalDescription"]["material"]);
        OpenMagnetics::CoreWrapper core(coreData);
        OpenMagnetics::CoilWrapper winding(windingData);
        OpenMagnetics::OperatingPoint operatingPoint(operatingPointData);
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");
        magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();
    }

    TEST(Test_Inductance_High_Flux_40_Web) {
        settings->reset();
        OpenMagnetics::clear_databases();

        // This tests checks that the operating is not crashing

        json coreData = json::parse(R"({"functionalDescription": {"gapping": [{"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 0.001, "sectionDimensions": null, "shape": null, "type": "subtractive"}, {"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}, {"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}], "material": "High Flux 40", "numberStacks": 1, "shape": {"aliases": [], "dimensions": {"A": 0.0391, "B": 0.0198, "C": 0.0125, "D": 0.0146, "E": 0.030100000000000002, "F": 0.0125, "G": 0.0, "H": 0.0}, "family": "etd", "familySubtype": "1", "magneticCircuit": null, "name": "ETD 39/20/13", "type": "standard"}, "type": "two-piece set"}, "geometricalDescription": null, "manufacturerInfo": null, "name": "My Core", "processedDescription": null})");
        json windingData =
            json::parse(R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary", "numberParallels": 1, "numberTurns": 24, "wire": "Dummy"}], "layersDescription": null, "sectionsDescription": null, "turnsDescription": null})");
        json inputsData = json::parse(R"({"designRequirements": {"altitude": null, "cti": null, "insulationType": null, "leakageInductance": null, "magnetizingInductance": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.0001279222825940401}, "name": null, "operatingTemperature": null, "overvoltageCategory": null, "pollutionDegree": null, "turnsRatios": []}, "operatingPoints": [{"conditions": {"ambientRelativeHumidity": null, "ambientTemperature": 25.0, "cooling": null, "name": null}, "excitationsPerWinding": [{"current": {"harmonics": null, "processed": null, "waveform": {"ancillaryLabel": null, "data": [-5.0, 5.0, -5.0], "numberPeriods": null, "time": [0.0, 2.5e-06, 1e-05]}}, "frequency": 100000.0, "magneticFieldStrength": null, "magneticFluxDensity": null, "magnetizingCurrent": null, "name": "My Operating Point", "voltage": null}], "name": null}]})");

        OpenMagnetics::CoreWrapper core(coreData);
        OpenMagnetics::CoilWrapper winding(windingData);
        OpenMagnetics::InputsWrapper inputs(inputsData);
        auto operatingPoint = inputs.get_operating_point(0);
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");
        magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();

    }

    TEST(Test_Inductance_Ferrite_Spacer) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 666;
        double frequency = 20000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_spacer_gap(0.003);
        double expectedValue = 13.5e-3;

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::CoilWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, -1, gapping, coreShape,
                                coreMaterial, core, winding, inputs);

        auto operatingPoint = inputs.get_operating_point(0);
        auto aux = magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint);
        double magnetizingInductance = aux.get_magnetizing_inductance().get_nominal().value();

        CHECK_CLOSE(expectedValue, magnetizingInductance, max_error * expectedValue);
    }

    TEST(Test_Inductance_Ferrite_Ground_Few_Turns) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 42;
        double numberTurns = 9;
        double frequency = 100000;
        std::string coreShape = "E 47/20/16";
        std::string coreMaterial = "N87";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0004);
        double expectedValue = 63e-6;

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::CoilWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, -1, gapping, coreShape,
                                coreMaterial, core, winding, inputs);

        auto operatingPoint = inputs.get_operating_point(0);
        double magnetizingInductance =
            magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();

        CHECK_CLOSE(expectedValue, magnetizingInductance, max_error * expectedValue);
    }

    TEST(Test_Inductance_Powder) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double dcCurrent = 96;
        double ambientTemperature = 25;
        double numberTurns = 13;
        double frequency = 68000;
        std::string coreShape = "E 42/21/15";
        std::string coreMaterial = "Edge 60";
        auto gapping = OpenMagneticsTesting::get_residual_gap();
        double expectedValue = 15.7e-6;

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::CoilWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, 20e6, gapping, coreShape,
                                coreMaterial, core, winding, inputs);

        auto operatingPoint = inputs.get_operating_point(0);
        auto aux = magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint);
        double magnetizingInductance = aux.get_magnetizing_inductance().get_nominal().value();

        CHECK_CLOSE(expectedValue, magnetizingInductance, max_error * expectedValue);
    }

    TEST(Test_NumberTurns_Ferrite_Ground) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 23.3e-3;
        double frequency = 20000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.003);

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::CoilWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedValue = 666;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, -1, desiredMagnetizingInductance, gapping,
                                coreShape, coreMaterial, core, winding, inputs);

        double numberTurns = magnetizing_inductance.calculate_number_turns_from_gapping_and_inductance(core, &inputs);

        CHECK_CLOSE(expectedValue, numberTurns, max_error * expectedValue);
    }

    TEST(Test_NumberTurns_Powder) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double dcCurrent = 96;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 15.7e-6;
        double frequency = 68000;
        std::string coreShape = "E 42/21/15";
        std::string coreMaterial = "Edge 60";
        auto gapping = OpenMagneticsTesting::get_residual_gap();

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::CoilWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedValue = 13;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, -1, desiredMagnetizingInductance, gapping,
                                coreShape, coreMaterial, core, winding, inputs);

        double numberTurns = magnetizing_inductance.calculate_number_turns_from_gapping_and_inductance(core, &inputs);

        CHECK_CLOSE(expectedValue, numberTurns, max_error * expectedValue);
    }

    TEST(Test_Gapping_Ferrite_Ground) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 23.3e-3;
        double numberTurns = 666;
        double frequency = 20000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::CoilWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedValue = 0.003;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, desiredMagnetizingInductance, {},
                                coreShape, coreMaterial, core, winding, inputs);

        auto gapping = magnetizing_inductance.calculate_gapping_from_number_turns_and_inductance(
            core, winding, &inputs, OpenMagnetics::GappingType::GROUND);

        CHECK_CLOSE(expectedValue, gapping[0].get_length(), max_error * expectedValue);
    }

    TEST(Test_Gapping_U_Shape_Ferrite_Ground) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 23.3e-3;
        double numberTurns = 666;
        double frequency = 20000;
        std::string coreShape = "U 26/22/16";
        std::string coreMaterial = "3C97";

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::CoilWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedValue = 0.0066;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, desiredMagnetizingInductance, {},
                                coreShape, coreMaterial, core, winding, inputs);

        auto gapping = magnetizing_inductance.calculate_gapping_from_number_turns_and_inductance(
            core, winding, &inputs, OpenMagnetics::GappingType::GROUND);

        CHECK_CLOSE(expectedValue, gapping[0].get_length(), max_error * expectedValue);
    }

    TEST(Test_Gapping_Ferrite_Distributed) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 23.3e-3;
        double numberTurns = 666;
        double frequency = 20000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::CoilWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedValue = 0.0004;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, desiredMagnetizingInductance, {},
                                coreShape, coreMaterial, core, winding, inputs);

        auto gapping = magnetizing_inductance.calculate_gapping_from_number_turns_and_inductance(
            core, winding, &inputs, OpenMagnetics::GappingType::DISTRIBUTED);

        CHECK_CLOSE(expectedValue, gapping[0].get_length(), max_error * expectedValue);
        CHECK_EQUAL(7UL, gapping.size());
    }

    TEST(Test_Gapping_Ferrite_Distributed_More_Gap_Precision) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 23.3e-3;
        double numberTurns = 666;
        double frequency = 20000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::CoilWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedValue = 0.00039;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, desiredMagnetizingInductance, {},
                                coreShape, coreMaterial, core, winding, inputs);

        auto gapping = magnetizing_inductance.calculate_gapping_from_number_turns_and_inductance(
            core, winding, &inputs, OpenMagnetics::GappingType::DISTRIBUTED, 5);

        CHECK_EQUAL(expectedValue, gapping[0].get_length());
        CHECK_EQUAL(7UL, gapping.size());
    }

    TEST(Test_Gapping_Classic_Web) {
        settings->reset();
        OpenMagnetics::clear_databases();

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
        OpenMagnetics::GappingType gappingType =
            magic_enum::enum_cast<OpenMagnetics::GappingType>("DISTRIBUTED").value();

        OpenMagnetics::CoreWrapper core(coreData);
        OpenMagnetics::CoilWrapper winding(windingData);
        OpenMagnetics::InputsWrapper inputs(inputsData);
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("CLASSIC");
        auto gapping =
            magnetizing_inductance.calculate_gapping_from_number_turns_and_inductance(core, winding, &inputs, gappingType, 5);

        CHECK(gapping.size() == 5);
    }

    TEST(Test_Gapping_Web) {
        settings->reset();
        OpenMagnetics::clear_databases();

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
        OpenMagnetics::GappingType gappingType = magic_enum::enum_cast<OpenMagnetics::GappingType>("GROUND").value();

        OpenMagnetics::CoreWrapper core(coreData);
        OpenMagnetics::CoilWrapper winding(windingData);
        OpenMagnetics::InputsWrapper inputs(inputsData);
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("CLASSIC");
        auto gapping =
            magnetizing_inductance.calculate_gapping_from_number_turns_and_inductance(core, winding, &inputs, gappingType, 5);
    }

    TEST(Test_Magnetizing_Inductance) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 42;
        double frequency = 20000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::CoilWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedInductanceValue = 215e-6;
        double currentPeakToPeak = 20;
        double voltagePeakToPeak = 2 * std::numbers::pi * frequency * expectedInductanceValue * currentPeakToPeak;

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
        OpenMagnetics::OperatingPointExcitation primaryExcitation =
            OpenMagnetics::InputsWrapper::get_primary_excitation(operatingPoint);

        CHECK_CLOSE(expectedInductanceValue, magnetizingInductance, max_error * expectedInductanceValue);
        CHECK_CLOSE(expectedMagneticFluxDensity, magneticFluxDensityWaveformPeak,
                    max_error * expectedMagneticFluxDensity);
        CHECK(bool(primaryExcitation.get_voltage()));
        CHECK(bool(primaryExcitation.get_magnetizing_current()));

        if (primaryExcitation.get_current()) {
            auto currentProcessed = primaryExcitation.get_current().value().get_processed().value();
            auto magnetizingCurrentProcessed = primaryExcitation.get_current().value().get_processed().value();
            CHECK_CLOSE(currentPeakToPeak,
                        operatingPoint.get_mutable_excitations_per_winding()[0]
                            .get_magnetizing_current()
                            .value()
                            .get_processed()
                            .value()
                            .get_peak_to_peak()
                            .value(),
                        max_error * currentPeakToPeak);
        }
    }

    TEST(Test_Gapping_Web_No_Voltage) {
        settings->reset();
        OpenMagnetics::clear_databases();

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
        OpenMagnetics::GappingType gappingType = magic_enum::enum_cast<OpenMagnetics::GappingType>("GROUND").value();

        OpenMagnetics::CoreWrapper core(coreData);
        OpenMagnetics::CoilWrapper winding(windingData);
        OpenMagnetics::InputsWrapper inputs(inputsData);
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("CLASSIC");
        auto gapping =
            magnetizing_inductance.calculate_gapping_from_number_turns_and_inductance(core, winding, &inputs, gappingType, 5);
        auto primaryExcitation = inputs.get_operating_point(0).get_mutable_excitations_per_winding()[0];
        double currentPeakToPeak = 10;

        CHECK(bool(primaryExcitation.get_voltage()));
        CHECK(bool(primaryExcitation.get_current()));
        CHECK(bool(primaryExcitation.get_magnetizing_current()));

        if (primaryExcitation.get_current()) {
            auto currentProcessed = primaryExcitation.get_current().value().get_processed().value();
            auto magnetizingCurrentProcessed = primaryExcitation.get_current().value().get_processed().value();
            CHECK_CLOSE(currentPeakToPeak,
                        inputs.get_operating_point(0)
                            .get_mutable_excitations_per_winding()[0]
                            .get_magnetizing_current()
                            .value()
                            .get_processed()
                            .value()
                            .get_peak_to_peak()
                            .value(),
                        max_error * currentPeakToPeak);
            CHECK_CLOSE(currentPeakToPeak,
                        inputs.get_operating_point(0)
                            .get_mutable_excitations_per_winding()[0]
                            .get_current()
                            .value()
                            .get_processed()
                            .value()
                            .get_peak_to_peak()
                            .value(),
                        max_error * currentPeakToPeak);
        }
    }

    TEST(Test_Inductance_Ferrite_Web_No_Voltage) {
        settings->reset();
        OpenMagnetics::clear_databases();

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

        OpenMagnetics::CoreWrapper core(coreData);
        OpenMagnetics::CoilWrapper winding(windingData);
        OpenMagnetics::OperatingPoint operatingPoint(operatingPointData);
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");
        magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();
        auto primaryExcitation = operatingPoint.get_mutable_excitations_per_winding()[0];
        double currentPeakToPeak = 10;
        double voltagePeakToPeak = 105;

        CHECK(bool(primaryExcitation.get_voltage()));
        CHECK(bool(primaryExcitation.get_current()));
        CHECK(bool(primaryExcitation.get_magnetizing_current()));

        if (primaryExcitation.get_current()) {
            auto currentProcessed = primaryExcitation.get_current().value().get_processed().value();
            auto magnetizingCurrentProcessed = primaryExcitation.get_current().value().get_processed().value();
            CHECK_CLOSE(currentPeakToPeak,
                        operatingPoint.get_mutable_excitations_per_winding()[0]
                            .get_magnetizing_current()
                            .value()
                            .get_processed()
                            .value()
                            .get_peak_to_peak()
                            .value(),
                        max_error * currentPeakToPeak);
            CHECK_CLOSE(currentPeakToPeak,
                        operatingPoint.get_mutable_excitations_per_winding()[0]
                            .get_current()
                            .value()
                            .get_processed()
                            .value()
                            .get_peak_to_peak()
                            .value(),
                        max_error * currentPeakToPeak);
        }
        if (primaryExcitation.get_voltage()) {
            CHECK_CLOSE(voltagePeakToPeak,
                        operatingPoint.get_mutable_excitations_per_winding()[0]
                            .get_voltage()
                            .value()
                            .get_processed()
                            .value()
                            .get_peak_to_peak()
                            .value(),
                        max_error * voltagePeakToPeak);
        }
    }

    TEST(Test_Magnetizing_Inductance_Toroid) {
        settings->reset();
        OpenMagnetics::clear_databases();

        settings->reset();
        OpenMagnetics::clear_databases();
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 42;
        double frequency = 20000;
        std::string coreShape = "T 58/41/18";
        std::string coreMaterial = "3C95";
        std::vector<OpenMagnetics::CoreGap> gapping = {};

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::CoilWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedValue = 6.6e-3;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, -1, gapping, coreShape,
                                coreMaterial, core, winding, inputs);

        auto operatingPoint = inputs.get_operating_point(0);
        double magnetizingInductance =
            magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();

        CHECK_CLOSE(expectedValue, magnetizingInductance, max_error * expectedValue);
    }

    TEST(Test_Magnetizing_Inductance_Toroid_Stacks) {
        settings->reset();
        OpenMagnetics::clear_databases();

        settings->reset();
        OpenMagnetics::clear_databases();
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 42;
        double frequency = 20000;
        std::string coreShape = "T 58/41/18";
        std::string coreMaterial = "3C95";
        std::vector<OpenMagnetics::CoreGap> gapping = {};

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::CoilWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedValue = 6.6e-3;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, -1, gapping, coreShape,
                                coreMaterial, core, winding, inputs);

        auto operatingPoint = inputs.get_operating_point(0);
        double magnetizingInductance =
            magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();

        CHECK_CLOSE(expectedValue, magnetizingInductance, max_error * expectedValue);

        core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, 2, coreMaterial);
        double magnetizingInductance2Stacks =
            magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();

        expectedValue = magnetizingInductance * 2;

        CHECK_CLOSE(expectedValue, magnetizingInductance2Stacks, max_error * expectedValue);
    }

    TEST(Test_Magnetizing_Inductance_RM14_20) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 29;
        double frequency = 100000;
        std::string coreShape = "RM 14/20";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);

        OpenMagnetics::CoreWrapper core = json::parse(R"({"name": "My Core", "functionalDescription": {"type": "two-piece set", "material": "3C97", "shape": {"aliases": ["RM 14LP", "RM 14/ILP", "RM 14/LP"], "dimensions": {"A": {"minimum": 0.0408, "maximum": 0.0422 }, "B": {"minimum": 0.010150000000000001, "maximum": 0.01025 }, "C": {"minimum": 0.018400000000000003, "maximum": 0.019000000000000003 }, "D": {"minimum": 0.00555, "maximum": 0.00585 }, "E": {"minimum": 0.029, "maximum": 0.0302 }, "F": {"minimum": 0.014400000000000001, "maximum": 0.015000000000000001 }, "G": {"minimum": 0.017 }, "H": {"minimum": 0.0054, "maximum": 0.005600000000000001 }, "J": {"minimum": 0.0335, "maximum": 0.0347 }, "R": {"maximum": 0.00030000000000000003 } }, "family": "rm", "familySubtype": "3", "name": "RM 14/20", "type": "standard", "magneticCircuit": "open"}, "gapping": [{"type": "subtractive", "length": 0.001 }, {"length": 0.000005, "type": "residual"}, {"length": 0.000005, "type": "residual"}, {"length": 0.000005, "type": "residual"} ], "numberStacks": 1 }, "geometricalDescription": null, "processedDescription": null })");
        OpenMagnetics::CoilWrapper winding; 
        OpenMagnetics::InputsWrapper inputs; 
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, -1, gapping, coreShape,
                                coreMaterial, core, winding, inputs);

        auto operatingPoint = inputs.get_operating_point(0);
        magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();
    }

    TEST(Test_Magnetizing_Inductance_Error_Web_0) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 10;
        double frequency = 20000;
        std::string coreShape = "E 65/32/27";
        std::string coreMaterial = "N95";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::CoilWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");

        double expectedValue = 19e-6;

        int numberStacks = 2;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, -1, gapping, coreShape,
                                coreMaterial, core, winding, inputs, 20, numberStacks);

        auto operatingPoint = inputs.get_operating_point(0);
        double magnetizingInductance =
            magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();

        CHECK_CLOSE(expectedValue, magnetizingInductance, max_error * expectedValue);
    }

    TEST(Test_Inductance_Powder_E_65) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double max_error = 0.15;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 10;
        double frequency = 100000;
        std::string coreShape = "E 65/32/27";
        std::string coreMaterial = "Kool Mµ 40";
        auto gapping = OpenMagneticsTesting::get_residual_gap();
        double expectedValue = 23e-6;

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::CoilWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, 20e6, gapping, coreShape,
                                coreMaterial, core, winding, inputs);

        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->set_magnetizing_inductance_include_air_inductance(true);

        auto operatingPoint = inputs.get_operating_point(0);
        auto aux = magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint);
        double magnetizingInductance = aux.get_magnetizing_inductance().get_nominal().value();

        CHECK_CLOSE(expectedValue, magnetizingInductance, max_error * expectedValue);
        settings->reset();
    }


    TEST(Test_Inductance_Powder_E_34) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double numberTurns = 10;
        double frequency = 100000;
        std::string coreShape = "E 34/14/9";
        std::string coreMaterial = "Edge 26";
        auto gapping = OpenMagneticsTesting::get_residual_gap();
        double expectedValue = 5.6e-6;

        OpenMagnetics::CoreWrapper core;
        OpenMagnetics::CoilWrapper winding;
        OpenMagnetics::InputsWrapper inputs;
        OpenMagnetics::MagnetizingInductance magnetizing_inductance("ZHANG");

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, 20e6, gapping, coreShape,
                                coreMaterial, core, winding, inputs);

        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->set_magnetizing_inductance_include_air_inductance(true);

        auto operatingPoint = inputs.get_operating_point(0);
        auto aux = magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint);
        double magnetizingInductance = aux.get_magnetizing_inductance().get_nominal().value();

        CHECK_CLOSE(expectedValue, magnetizingInductance, max_error * expectedValue);
        settings->reset();
    }


}