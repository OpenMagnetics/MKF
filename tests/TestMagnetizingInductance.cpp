#include <source_location>
#include "processors/CircuitSimulatorInterface.h"
#include "processors/ResultsAnalyzer.h"
#include "processors/MagneticSimulator.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/InitialPermeability.h"
#include "constructive_models/CorePiece.h"
#include "physical_models/Reluctance.h"
#include "constructive_models/Bobbin.h"
#include "constructive_models/Magnetic.h"
#include "TestingUtils.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <chrono>
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
            "magneticCircuit": null, "name": "ETD 39/20/13", "type": "standard"}, "type": "twoPieceSet"}})");
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
        double computedMagnetizingInductance = magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();
        // The computed magnetizing inductance must be finite and positive.
        CHECK(std::isfinite(computedMagnetizingInductance));
        CHECK(computedMagnetizingInductance > 0);
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
            "type": "standard"}, "type": "twoPieceSet"}, "geometricalDescription": null,
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
        double computedMagnetizingInductance = magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();
        // The computed magnetizing inductance must be finite and positive.
        CHECK(std::isfinite(computedMagnetizingInductance));
        CHECK(computedMagnetizingInductance > 0);
    }

    TEST_CASE("Test_Inductance_High_Flux_40_Web", "[physical-model][magnetizing-inductance][smoke-test]") {
        settings.reset();
        clear_databases();

        // This tests checks that the operating is not crashing

        json coreData = json::parse(R"({"functionalDescription": {"gapping": [{"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 0.001, "sectionDimensions": null, "shape": null, "type": "subtractive"}, {"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}, {"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}], "material": "High Flux 40", "numberStacks": 1, "shape": {"aliases": [], "dimensions": {"A": 0.0391, "B": 0.0198, "C": 0.0125, "D": 0.0146, "E": 0.030100000000000002, "F": 0.0125, "G": 0.0, "H": 0.0}, "family": "etd", "familySubtype": "1", "magneticCircuit": null, "name": "ETD 39/20/13", "type": "standard"}, "type": "twoPieceSet"}, "geometricalDescription": null, "manufacturerInfo": null, "name": "My Core", "processedDescription": null})");
        json windingData =
            json::parse(R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary", "numberParallels": 1, "numberTurns": 24, "wire": "Dummy"}], "layersDescription": null, "sectionsDescription": null, "turnsDescription": null})");
        json inputsData = json::parse(R"({"designRequirements": {"altitude": null, "cti": null, "insulationType": null, "leakageInductance": null, "magnetizingInductance": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.0001279222825940401}, "name": null, "operatingTemperature": null, "overvoltageCategory": null, "pollutionDegree": null, "turnsRatios": []}, "operatingPoints": [{"conditions": {"ambientRelativeHumidity": null, "ambientTemperature": 25.0, "cooling": null, "name": null}, "excitationsPerWinding": [{"current": {"harmonics": null, "processed": null, "waveform": {"ancillaryLabel": null, "data": [-5.0, 5.0, -5.0], "numberPeriods": null, "time": [0.0, 2.5e-06, 1e-05]}}, "frequency": 100000.0, "magneticFieldStrength": null, "magneticFluxDensity": null, "magnetizingCurrent": null, "name": "My Operating Point", "voltage": null}], "name": null}]})");

        Core core(coreData);
        OpenMagnetics::Coil winding(windingData);
        OpenMagnetics::Inputs inputs(inputsData);
        auto operatingPoint = inputs.get_operating_point(0);
        MagnetizingInductance magnetizing_inductance("ZHANG");
        double computedMagnetizingInductance = magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();
        // The computed magnetizing inductance must be finite and positive.
        CHECK(std::isfinite(computedMagnetizingInductance));
        CHECK(computedMagnetizingInductance > 0);
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

        double numberTurns = magnetizing_inductance.calculate_number_turns_from_gapping_and_inductance(core, winding, &inputs);

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

        double numberTurns = magnetizing_inductance.calculate_number_turns_from_gapping_and_inductance(core, winding, &inputs);

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(numberTurns, max_error * expectedValue));
    }

    TEST_CASE("Test_NumberTurns_Nominal_Picks_Closest_Not_Ceil", "[physical-model][magnetizing-inductance][bug]") {
        // ABT #600: with a NOMINAL target the recommender ceil()-ed against the
        // target, accepting a +26.6% overshoot to avoid a -0.6% undershoot
        // (E 42/21/15 / 3C95 / 1 mm gap / 10 uH: L(8)=9.94 uH, L(9)=12.58 uH,
        // and it returned 9). NOMINAL must pick the neighbour with the smallest
        // absolute error; MINIMUM keeps floor-clearing semantics.
        settings.reset();
        clear_databases();

        double dcCurrent = 0;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 10e-6;
        double frequency = 100000;
        std::string coreShape = "E 42/21/15";
        std::string coreMaterial = "3C95";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);

        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, -1, desiredMagnetizingInductance, gapping,
                                coreShape, coreMaterial, core, winding, inputs);

        int numberTurns = magnetizing_inductance.calculate_number_turns_from_gapping_and_inductance(core, winding, &inputs);

        auto inductanceAt = [&](int n) {
            winding.get_mutable_functional_description()[0].set_number_turns(n);
            auto operatingPoint = inputs.get_operating_point(0);
            return resolve_dimensional_values(
                magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint)
                    .get_magnetizing_inductance());
        };

        // Neither neighbour may deliver an inductance closer to the target.
        double errorAtN = std::fabs(inductanceAt(numberTurns) - desiredMagnetizingInductance);
        REQUIRE(errorAtN <= std::fabs(inductanceAt(numberTurns + 1) - desiredMagnetizingInductance));
        if (numberTurns > 1) {
            REQUIRE(errorAtN <= std::fabs(inductanceAt(numberTurns - 1) - desiredMagnetizingInductance));
        }

        // The same target requested as a hard floor (MINIMUM) must clear it.
        int numberTurnsMinimum = magnetizing_inductance.calculate_number_turns_from_gapping_and_inductance(
            core, winding, &inputs, DimensionalValues::MINIMUM);
        REQUIRE(numberTurnsMinimum >= numberTurns);
        REQUIRE(inductanceAt(numberTurnsMinimum) >= desiredMagnetizingInductance * 0.999);
    }

    TEST_CASE("Test_NumberTurns_Legacy_NoCoil_Overload_Matches_SingleWinding",
              "[physical-model][magnetizing-inductance][smoke-test]") {
        // The legacy 3-argument overload (no coil) must reproduce the coil-aware
        // result for a single primary winding — it synthesizes exactly that coil
        // internally. Magnetizing inductance is N^2 / R_core (primary-referred,
        // wire-independent), so the two paths must agree to integer rounding.
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

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, -1, desiredMagnetizingInductance, gapping,
                                coreShape, coreMaterial, core, winding, inputs);

        int withCoil = magnetizing_inductance.calculate_number_turns_from_gapping_and_inductance(core, winding, &inputs);
        int legacyNoCoil = magnetizing_inductance.calculate_number_turns_from_gapping_and_inductance(core, &inputs);

        REQUIRE(legacyNoCoil == withCoil);
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

    // ABT #1093: the DC-bias permeability used to come from a fixed-point iteration that
    // diverges at the B-H knee; with 2 A of bias on a 15-turn P 11/7/I the search ran away to
    // a saturated permeability and returned the residual gap for any target above ~20 µH.
    // The bias is now solved on the material's own B-H curve, so the gap is real and the
    // bias-aware inductance of the gapped core lands on the target.
    TEST_CASE("Test_Gapping_With_Dc_Bias_Does_Not_Collapse_To_Residual", "[physical-model][magnetizing-inductance][gapping][dc-bias][smoke-test]") {
        settings.reset();
        clear_databases();
        double dcCurrent = 2;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 36.5e-6;
        double numberTurns = 15;
        double frequency = 100000;
        std::string coreShape = "P 11/7/I";
        std::string coreMaterial = "98";
        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");
        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, desiredMagnetizingInductance, {},
                                coreShape, coreMaterial, core, winding, inputs, 1);

        auto gapping = magnetizing_inductance.calculate_gapping_from_number_turns_and_inductance(
            core, winding, &inputs, GappingType::GROUND, 6);
        UNSCOPED_INFO("ground gap " << gapping[0].get_length() * 1e6 << " um");
        REQUIRE(gapping[0].get_length() > 50e-6);   // not the 5 um residual
        REQUIRE(gapping[0].get_length() < 400e-6);

        Core gappedCore(core);
        gappedCore.get_mutable_functional_description().set_gapping(gapping);
        gappedCore.set_processed_description(std::nullopt);
        gappedCore.process_data();
        gappedCore.process_gap();
        auto operatingPoint = inputs.get_operating_point(0);
        double inductance = magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(gappedCore, winding, &operatingPoint)
                                .get_magnetizing_inductance().get_nominal().value();
        UNSCOPED_INFO("bias-aware inductance " << inductance * 1e6 << " uH for target " << desiredMagnetizingInductance * 1e6);
        REQUIRE_THAT(inductance, Catch::Matchers::WithinRel(desiredMagnetizingInductance, 0.25));
        settings.reset();
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

        // Re-pinned 2026-07-30 (ABT #378, user-approved): 0.0004 -> 0.00039 after Zhang's h was
        // read as the paper defines it (Fig. 7: "2h is the height of a segment of core limb"),
        // dropping the old clamp of h up to the column WIDTH. Less modelled fringing means more
        // reluctance per unit of gap, so a slightly SHORTER gap now reaches the same 23.3 mH —
        // the shift is in the physically expected direction, and at 1e-5 m it is ten times the
        // solver's 1e-6 rounding quantum, so it is a real change rather than noise. This is a
        // solved SOLVER OUTPUT, not measured data.
        // NOTE: this fixture is DISTRIBUTED gapping, which is exactly the case where MKF's h is
        // still over-estimated — it measures to the limb end, while the paper's h is half the
        // core segment SHARED with the neighbouring gap. Expect this value to move once that
        // refinement lands (tracked on ABT #378).
        double expectedValue = 0.00039;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, desiredMagnetizingInductance, {},
                                coreShape, coreMaterial, core, winding, inputs);

        auto gapping = magnetizing_inductance.calculate_gapping_from_number_turns_and_inductance(
            core, winding, &inputs, GappingType::DISTRIBUTED, 5);

        REQUIRE(expectedValue == gapping[0].get_length());
        REQUIRE(7UL == gapping.size());
    }

    // ABT #1093: these three web-shaped cases used to carry a 46 A DC current (41..51 A) that put
    // 1.7 T, 765 T and 7.6 T in 3C95 (saturation 0.53 T); the old fixed-point bias iteration
    // silently produced a saturated permeability and some gap. The bias is now solved on the
    // material's magnetisation curve and a flux density above saturation is refused (see
    // Test_Gapping_Bias_Beyond_Saturation_Throws), so the currents are 2 A DC with the same
    // 10 A peak-to-peak, and Test_Gapping_Web asks 4.65 µH of its single turn instead of 4.65 mH
    // (a thousand times the ungapped AL of an ETD 54). The cases still exercise what they were
    // written for: the web JSON goes through distributed / ground gapping without crashing.
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
            "open", "name": "ETD 54/28/19", "type": "standard"}, "type": "twoPieceSet"},
            "geometricalDescription": [{"coordinates": [0.0, 0.0, 0.0], "dimensions": null,
            "insulationMaterial": null, "machining": [{"coordinates": [0.0, 7.5e-05, 0.0], "length":
            0.00015}, {"coordinates": [0.0, 0.0073, 0.0], "length": 0.0003}], "material": "3C97",
            "rotation": [3.141592653589793, 3.141592653589793, 0.0], "shape": {"aliases": [], "dimensions":
            {"A": 0.0391, "B": 0.0198, "C": 0.0125, "D": 0.0146, "E": 0.030100000000000002, "F": 0.0125},
            "family": "etd", "familySubtype": "1", "magneticCircuit": null, "name": "ETD 39/20/13",
            "type": "standard"}, "type": "halfSet"}, {"coordinates": [0.0, 0.0, 0.0], "dimensions":
            null, "insulationMaterial": null, "machining": [{"coordinates": [0.0, -0.0073, 0.0], "length":
            0.0003}, {"coordinates": [0.0, -7.5e-05, 0.0], "length": 0.00015}], "material": "3C97",
            "rotation": [0.0, 0.0, 0.0], "shape": {"aliases": [], "dimensions": {"A": 0.0391, "B":
            0.0198, "C": 0.0125, "D": 0.0146, "E": 0.030100000000000002, "F": 0.0125}, "family": "etd",
            "familySubtype": "1", "magneticCircuit": null, "name": "ETD 39/20/13", "type": "standard"},
            "type": "halfSet"}], "processedDescription": {"columns": [{"area": 0.000123, "coordinates":
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
            "data": [-3.0, 7.0, -3.0], "numberPeriods": null, "time": [0.0, 2.4999999999999998e-06, 1e-05]}},
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

    TEST_CASE("Test_Gapping_Bias_Beyond_Saturation_Throws", "[physical-model][magnetizing-inductance][gapping][dc-bias][smoke-test]") {
        settings.reset();
        clear_databases();
        // 46 A of DC through 40 turns asking 412 µH on an ETD 29 puts ~6 T in the ferrite: no gap
        // holds that inductance at that bias. Refused loudly instead of returning a gap that does
        // not deliver the target (ABT #1093).
        double dcCurrent = 46;
        double ambientTemperature = 25;
        double desiredMagnetizingInductance = 412.7e-6;
        double numberTurns = 40;
        double frequency = 100000;
        std::string coreShape = "ETD 29";
        std::string coreMaterial = "3C97";
        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");
        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, desiredMagnetizingInductance, {},
                                coreShape, coreMaterial, core, winding, inputs, 10);
        REQUIRE_THROWS_WITH(magnetizing_inductance.calculate_gapping_from_number_turns_and_inductance(core, winding, &inputs, GappingType::GROUND, 5),
                            Catch::Matchers::ContainsSubstring("above its saturation"));
        settings.reset();
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
            "open", "name": "ETD 54/28/19", "type": "standard"}, "type": "twoPieceSet"},
            "geometricalDescription": null, "processedDescription": null})");
        json windingData =
            json::parse(R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary",
                        "numberParallels": 1, "numberTurns": 1, "wire": "Dummy"}], "layersDescription":
                        null, "sectionsDescription": null, "turnsDescription": null})");
        json inputsData = json::parse(
            R"({"designRequirements": {"altitude": null, "cti": null, "insulationType": null,
            "leakageInductance": null, "magnetizingInductance": {"excludeMaximum": null, "excludeMinimum":
            null, "maximum": null, "minimum": null, "nominal": 4.654652816558039e-06}, "name": null,
            "operatingTemperature": null, "overvoltageCategory": null, "pollutionDegree": null,
            "turnsRatios": []}, "operatingPoints": [{"conditions": {"ambientRelativeHumidity": null,
            "ambientTemperature": 25.0, "cooling": null, "name": null}, "excitationsPerWinding":
            [{"current": {"harmonics": null, "processed": null, "waveform": {"ancillaryLabel": null,
            "data": [-3.0, 7.0, -3.0], "numberPeriods": null, "time": [0.0, 2.4999999999999998e-06, 1e-05]}},
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
        // Repro point: a ground gapping must be computable for these web inputs
        // (the function's 5th argument is rounding decimals, not a gap count; GROUND
        // gapping on this ETD yields one subtractive gap plus residual gaps).
        REQUIRE(gapping.size() > 0);
        for (auto& gap : gapping) {
            CHECK(gap.get_length() > 0);
        }
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
        // 10 A peak-to-peak puts B_peak = L I_pk / (N A_e) = 0.335 T in the 3C97, against B_sat ~ 0.53 T at
        // 25 C. It was 20 A, a 0.669 T peak the material cannot carry, which the inductance model now refuses
        // as a saturation (ABT #1740).
        double currentPeakToPeak = 10;
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

        // A web request: ETD 54 3C95 with ONE turn, 46.5 uH, current-only excitation -3..7 A at 100 kHz.
        // No gap gives 46.5 uH with one turn on an unsaturated core: B pk = 1.16 T against
        // B_sat 0.53 T (3C95, 25 C). The gap solve must REFUSE with CALCULATION_CORE_SATURATES naming
        // B and B_sat, never return an air-core gap as a design. (Formerly asserted a returned gap and
        // an induced voltage; rewritten to expect the refusal, Alf approved 2026-10-14, ABT #1855.)
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
            "open", "name": "ETD 54/28/19", "type": "standard"}, "type": "twoPieceSet"},
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
            "data": [-3.0, 7.0, -3.0], "numberPeriods": null, "time": [0.0, 2.5e-06, 1e-05]}}, "frequency":
            100000.0, "magneticField": null, "magneticFluxDensity": null, "magnetizingCurrent": null, "name":
            "My Operating Point"}], "name": null}]})");
        GappingType gappingType = magic_enum::enum_cast<GappingType>("GROUND").value();

        Core core(coreData);
        OpenMagnetics::Coil winding(windingData);
        OpenMagnetics::Inputs inputs(inputsData);
        MagnetizingInductance magnetizing_inductance("CLASSIC");
        bool refused = false;
        try {
            magnetizing_inductance.calculate_gapping_from_number_turns_and_inductance(core, winding, &inputs, gappingType, 5);
        }
        catch (const SaturationException& e) {
            refused = true;
            CHECK(e.code() == OpenMagnetics::ErrorCode::CALCULATION_CORE_SATURATES);
            const std::string message = e.what();
            CHECK_THAT(message, Catch::Matchers::ContainsSubstring("CALCULATION_CORE_SATURATES"));
            CHECK_THAT(message, Catch::Matchers::ContainsSubstring("B peak = 1.16"));
            CHECK_THAT(message, Catch::Matchers::ContainsSubstring("B_sat = 0.53"));
            CHECK_THAT(message, Catch::Matchers::ContainsSubstring("3C95"));
        }
        REQUIRE(refused);
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
            "magneticCircuit": null, "name": "ETD 39/20/13", "type": "standard"}, "type": "twoPieceSet"}})");
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
        double computedMagnetizingInductance = magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();
        // The computed magnetizing inductance must be finite and positive.
        CHECK(std::isfinite(computedMagnetizingInductance));
        CHECK(computedMagnetizingInductance > 0);
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

    TEST_CASE("Test_Magnetizing_Inductance_Gapless_Core_Has_No_Gap_Fields", "[physical-model][magnetizing-inductance][reluctance][bug][abt-1714]") {
        // ABT #1714: a gapless core (here a toroid) got gappingReluctance = 0, maximumFringingFactor = 1,
        // maximumStorableMagneticEnergyGapping = 0 and reluctancePerGap = [] in its magnetizing
        // inductance output. MAS gappingReluctance has exclusiveMinimum 0, so every exported MAS of a
        // toroid failed schema validation. With no gap those fields describe nothing: they are absent.
        settings.reset();
        clear_databases();
        double numberTurns = 42;
        std::vector<CoreGap> noGaps = {};
        Core toroid;
        OpenMagnetics::Coil toroidWinding;
        OpenMagnetics::Inputs toroidInputs;
        prepare_test_parameters(0, 25, 20000, numberTurns, -1, noGaps, "T 58/41/18", "3C95", toroid, toroidWinding, toroidInputs);
        auto toroidOperatingPoint = toroidInputs.get_operating_point(0);
        MagnetizingInductance magnetizingInductanceModel("ZHANG");
        auto toroidOutput = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(toroid, toroidWinding, &toroidOperatingPoint);
        // MAS.hpp serialises an unset optional as null (as it does measurementCondition); the MAS
        // exporters drop nulls, so null here is an absent field in the exported file.
        json toroidJson;
        to_json(toroidJson, toroidOutput);
        INFO(toroidJson.dump());
        CHECK(toroidJson.value("gappingReluctance", json()).is_null());
        CHECK(toroidJson.value("maximumFringingFactor", json()).is_null());
        CHECK(toroidJson.value("maximumStorableMagneticEnergyGapping", json()).is_null());
        CHECK(toroidJson.value("reluctancePerGap", json()).is_null());
        REQUIRE(toroidJson.contains("coreReluctance"));
        CHECK(toroidJson.at("coreReluctance").get<double>() > 0);
        // No gap in series: the core reluctance is the material's alone.
        CHECK(toroidJson.at("coreReluctance").get<double>() == toroidJson.at("ungappedCoreReluctance").get<double>());
        CHECK(toroidOutput.get_magnetizing_inductance().get_nominal().value() > 0);

        // A gapped core keeps every gap field, and each is a positive number.
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        Core gapped;
        OpenMagnetics::Coil gappedWinding;
        OpenMagnetics::Inputs gappedInputs;
        prepare_test_parameters(0, 25, 20000, numberTurns, -1, gapping, "ETD 29", "3C95", gapped, gappedWinding, gappedInputs);
        auto gappedOperatingPoint = gappedInputs.get_operating_point(0);
        auto gappedOutput = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(gapped, gappedWinding, &gappedOperatingPoint);
        REQUIRE(gappedOutput.get_gapping_reluctance());
        CHECK(gappedOutput.get_gapping_reluctance().value() > 0);
        REQUIRE(gappedOutput.get_maximum_fringing_factor());
        CHECK(gappedOutput.get_maximum_fringing_factor().value() > 0);
        REQUIRE(gappedOutput.get_maximum_storable_magnetic_energy_gapping());
        CHECK(gappedOutput.get_maximum_storable_magnetic_energy_gapping().value() > 0);
        REQUIRE(gappedOutput.get_reluctance_per_gap());
        CHECK(!gappedOutput.get_reluctance_per_gap()->empty());
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

        Core core = json::parse(R"({"name": "My Core", "functionalDescription": {"type": "twoPieceSet", "material": "3C97", "shape": {"aliases": ["RM 14LP", "RM 14/ILP", "RM 14/LP"], "dimensions": {"A": {"minimum": 0.0408, "maximum": 0.0422 }, "B": {"minimum": 0.010150000000000001, "maximum": 0.01025 }, "C": {"minimum": 0.018400000000000003, "maximum": 0.019000000000000003 }, "D": {"minimum": 0.00555, "maximum": 0.00585 }, "E": {"minimum": 0.029, "maximum": 0.0302 }, "F": {"minimum": 0.014400000000000001, "maximum": 0.015000000000000001 }, "G": {"minimum": 0.017 }, "H": {"minimum": 0.0054, "maximum": 0.005600000000000001 }, "J": {"minimum": 0.0335, "maximum": 0.0347 }, "R": {"maximum": 0.00030000000000000003 } }, "family": "rm", "familySubtype": "3", "name": "RM 14/20", "type": "standard", "magneticCircuit": "open"}, "gapping": [{"type": "subtractive", "length": 0.001 }, {"length": 0.000005, "type": "residual"}, {"length": 0.000005, "type": "residual"}, {"length": 0.000005, "type": "residual"} ], "numberStacks": 1 }, "geometricalDescription": null, "processedDescription": null })");
        OpenMagnetics::Coil winding; 
        OpenMagnetics::Inputs inputs; 
        MagnetizingInductance magnetizing_inductance("ZHANG");

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, numberTurns, -1, gapping, coreShape,
                                coreMaterial, core, winding, inputs);

        auto operatingPoint = inputs.get_operating_point(0);
        double computedMagnetizingInductance = magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();
        // The computed magnetizing inductance must be finite and positive.
        CHECK(std::isfinite(computedMagnetizingInductance));
        CHECK(computedMagnetizingInductance > 0);
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

        Core core = json::parse(R"({"name": "650-4637", "functionalDescription": {"type": "twoPieceSet", "material": "TP4A", "shape": {"aliases": ["E 16/5", "EF 16"], "dimensions": {"A": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0167, "minimum": 0.0155, "nominal": null}, "B": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0082, "minimum": 0.0079, "nominal": null}, "C": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0047, "minimum": 0.0043, "nominal": null}, "D": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0061, "minimum": 0.0057, "nominal": null}, "E": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0119, "minimum": 0.0113, "nominal": null}, "F": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0047, "minimum": 0.0044, "nominal": null}}, "family": "e", "familySubtype": null, "magneticCircuit": "open", "name": "E 16/8/5", "type": "standard"}, "gapping": [], "numberStacks": 1}, "processedDescription": {"columns": [{"area": 2.1e-05, "coordinates": [0.0, 0.0, 0.0], "depth": 0.004501, "height": 0.011802, "minimumDepth": null, "minimumWidth": null, "shape": "rectangular", "type": "central", "width": 0.00455}, {"area": 1.1e-05, "coordinates": [0.006925, 0.0, 0.0], "depth": 0.004501, "height": 0.011802, "minimumDepth": null, "minimumWidth": null, "shape": "rectangular", "type": "lateral", "width": 0.002251}, {"area": 1.1e-05, "coordinates": [-0.006925, 0.0, 0.0], "depth": 0.004501, "height": 0.011802, "minimumDepth": null, "minimumWidth": null, "shape": "rectangular", "type": "lateral", "width": 0.002251}], "depth": 0.0045000000000000005, "effectiveParameters": {"effectiveArea": 2.0062091987236854e-05, "effectiveLength": 0.03756497447228765, "effectiveVolume": 7.53631973361239e-07, "minimumArea": 1.935000000000001e-05}, "height": 0.016100000000000003, "width": 0.0161, "windingWindows": [{"angle": null, "area": 4.1595e-05, "coordinates": [0.002275, 0.0], "height": 0.011800000000000001, "radialHeight": null, "sectionsAlignment": null, "sectionsOrientation": null, "shape": null, "width": 0.0035249999999999995}]}})");
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
        std::string coreMaterial = "Kool M\xC2\xB5 40";
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

        json coreData = json::parse(R"({"distributorsInfo":null,"functionalDescription":{"coating":null,"gapping":[],"material":"MPP 26","numberStacks":1,"shape":{"aliases":["R 80/20/50"],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.08},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.05},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.02}},"family":"t","familySubtype":null,"magneticCircuit":"closed","name":"T 80/20/50","type":"standard"},"type":"toroidal"},"geometricalDescription":[{"coordinates":[0.0,0.0,0.0],"dimensions":null,"insulationMaterial":null,"machining":null,"material":"MPP 26","rotation":[1.5707963267948966,1.5707963267948966,0.0],"shape":{"aliases":["R 80/20/50"],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.08},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.05},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.02}},"family":"t","familySubtype":null,"magneticCircuit":"closed","name":"T 80/20/50","type":"standard"},"type":"toroidal"}],"manufacturerInfo":null,"name":"custom","processedDescription":{"columns":[{"area":0.0003,"coordinates":[0.0,0.0,0.0],"depth":0.02,"height":0.20420352248333656,"minimumDepth":null,"minimumWidth":null,"shape":"rectangular","type":"central","width":0.015}],"depth":0.02,"effectiveParameters":{"effectiveArea":0.0003,"effectiveLength":0.20420352248333654,"effectiveVolume":6.126105674500096e-05,"minimumArea":0.0003},"height":0.08,"thermalResistance":null,"width":0.08,"windingWindows":[{"angle":360.0,"area":0.001963495408493621,"coordinates":[0.015,0.0],"height":null,"radialHeight":0.025,"sectionsAlignment":null,"sectionsOrientation":null,"shape":null,"width":null}]}})");
        json windingData = json::parse(R"({"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.01,"columnShape":"rectangular","columnThickness":0.0,"columnWidth":0.0075,"coordinates":[0.0,0.0,0.0],"pins":null,"wallThickness":0.0,"windingWindows":[{"angle":360.0,"area":0.001963495408493621,"coordinates":[0.025,0.0,0.0],"height":null,"radialHeight":0.025,"sectionsAlignment":"innerOrTop","sectionsOrientation":"overlapping","shape":"round","width":null}]}},"functionalDescription":[{"connections":null,"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":1,"wire":{"coating":{"breakdownVoltage":4600.0,"grade":2,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1.8834326265752323e-05},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.004897},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Nearson","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 4.5 - Heavy Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.004989},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"4.5 AWG","strand":null,"type":"round"}},{"connections":null,"isolationSide":"primary","name":"Secondary","numberParallels":1,"numberTurns":1,"wire":{"coating":{"breakdownVoltage":4600.0,"grade":2,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1.8834326265752323e-05},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.004897},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Nearson","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 4.5 - Heavy Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.004989},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"4.5 AWG","strand":null,"type":"round"}},{"connections":null,"isolationSide":"primary","name":"Tertiary","numberParallels":1,"numberTurns":1,"wire":{"coating":{"breakdownVoltage":4600.0,"grade":2,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1.8834326265752323e-05},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.004897},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Nearson","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 4.5 - Heavy Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.004989},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"4.5 AWG","strand":null,"type":"round"}}],"layersDescription":[{"additionalCoordinates":null,"coordinateSystem":"polar","coordinates":[0.0024945,180.0],"dimensions":[0.004989,360.0],"fillingFactor":0.03967937689697995,"insulationMaterial":null,"name":"Primary section 0 layer 0","orientation":"overlapping","partialWindings":[{"connections":null,"parallelsProportion":[1.0],"winding":"Primary"}],"section":"Primary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveParallels"},{"additionalCoordinates":[[-0.0200015,180.0]],"coordinateSystem":"polar","coordinates":[0.0050015,180.0],"dimensions":[2.5e-05,360.0],"fillingFactor":1.0,"insulationMaterial":null,"name":"Insulation between Primary and Primary section 1 layer 0","orientation":"overlapping","partialWindings":[],"section":"Insulation between Primary and Primary section 1","turnsAlignment":"spread","type":"insulation","windingStyle":null},{"additionalCoordinates":null,"coordinateSystem":"polar","coordinates":[0.007508500000000001,180.0],"dimensions":[0.004989,360.0],"fillingFactor":0.03967937689697995,"insulationMaterial":null,"name":"Secondary section 0 layer 0","orientation":"overlapping","partialWindings":[{"connections":null,"parallelsProportion":[1.0],"winding":"Secondary"}],"section":"Secondary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveParallels"},{"additionalCoordinates":[[-0.0250155,180.0]],"coordinateSystem":"polar","coordinates":[0.0100155,180.0],"dimensions":[2.5e-05,360.0],"fillingFactor":1.0,"insulationMaterial":null,"name":"Insulation between Secondary and Secondary section 3 layer 0","orientation":"overlapping","partialWindings":[],"section":"Insulation between Secondary and Secondary section 3","turnsAlignment":"spread","type":"insulation","windingStyle":null},{"additionalCoordinates":null,"coordinateSystem":"polar","coordinates":[0.0125225,180.0],"dimensions":[0.004989,360.0],"fillingFactor":0.03967937689697995,"insulationMaterial":null,"name":"Tertiary section 0 layer 0","orientation":"overlapping","partialWindings":[{"connections":null,"parallelsProportion":[1.0],"winding":"Tertiary"}],"section":"Tertiary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveParallels"},{"additionalCoordinates":[[-0.0300295,180.0]],"coordinateSystem":"polar","coordinates":[0.0150295,180.0],"dimensions":[2.5e-05,360.0],"fillingFactor":1.0,"insulationMaterial":null,"name":"Insulation between Tertiary and Tertiary section 5 layer 0","orientation":"overlapping","partialWindings":[],"section":"Insulation between Tertiary and Tertiary section 5","turnsAlignment":"spread","type":"insulation","windingStyle":null}],"sectionsDescription":[{"coordinateSystem":"polar","coordinates":[0.0024945,180.0],"dimensions":[0.004989,360.0],"fillingFactor":0.035281331722710724,"layersAlignment":null,"layersOrientation":"overlapping","margin":[0.0,0.0],"name":"Primary section 0","partialWindings":[{"connections":null,"parallelsProportion":[1.0],"winding":"Primary"}],"type":"conduction","windingStyle":"windByConsecutiveParallels"},{"coordinateSystem":"polar","coordinates":[0.005001500000000001,180.0],"dimensions":[2.5e-05,360.0],"fillingFactor":1.0,"layersAlignment":null,"layersOrientation":"overlapping","margin":null,"name":"Insulation between Primary and Primary section 1","partialWindings":[],"type":"insulation","windingStyle":null},{"coordinateSystem":"polar","coordinates":[0.007508500000000001,180.0],"dimensions":[0.004989,360.0],"fillingFactor":0.04539484956038453,"layersAlignment":null,"layersOrientation":"overlapping","margin":[0.0,0.0],"name":"Secondary section 0","partialWindings":[{"connections":null,"parallelsProportion":[1.0],"winding":"Secondary"}],"type":"conduction","windingStyle":"windByConsecutiveParallels"},{"coordinateSystem":"polar","coordinates":[0.010015500000000002,180.0],"dimensions":[2.5e-05,360.0],"fillingFactor":1.0,"layersAlignment":null,"layersOrientation":"overlapping","margin":null,"name":"Insulation between Secondary and Secondary section 3","partialWindings":[],"type":"insulation","windingStyle":null},{"coordinateSystem":"polar","coordinates":[0.0125225,180.0],"dimensions":[0.004989,360.0],"fillingFactor":0.06363646652658514,"layersAlignment":null,"layersOrientation":"overlapping","margin":[0.0,0.0],"name":"Tertiary section 0","partialWindings":[{"connections":null,"parallelsProportion":[1.0],"winding":"Tertiary"}],"type":"conduction","windingStyle":"windByConsecutiveParallels"},{"coordinateSystem":"polar","coordinates":[0.015029500000000003,180.0],"dimensions":[2.5e-05,360.0],"fillingFactor":1.0,"layersAlignment":null,"layersOrientation":"overlapping","margin":null,"name":"Insulation between Tertiary and Tertiary section 5","partialWindings":[],"type":"insulation","windingStyle":null}],"turnsDescription":[{"additionalCoordinates":[[-0.042494500000000004,5.204075340636721e-18]],"angle":null,"coordinateSystem":"cartesian","coordinates":[-0.0225055,2.756128853821076e-18],"dimensions":[0.004989,0.004989],"layer":"Primary section 0 layer 0","length":0.08567340574875948,"name":"Primary parallel 0 turn 0","orientation":"clockwise","parallel":0,"rotation":180.0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":[[-0.0475085,5.818113245729203e-18]],"angle":null,"coordinateSystem":"cartesian","coordinates":[-0.0174915,2.142090948728593e-18],"dimensions":[0.004989,0.004989],"layer":"Secondary section 0 layer 0","length":0.11717729687895795,"name":"Secondary parallel 0 turn 0","orientation":"clockwise","parallel":0,"rotation":180.0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":[[-0.0525225,6.432151150821686e-18]],"angle":null,"coordinateSystem":"cartesian","coordinates":[-0.0124775,1.52805304363611e-18],"dimensions":[0.004989,0.004989],"layer":"Tertiary section 0 layer 0","length":0.14868118800915636,"name":"Tertiary parallel 0 turn 0","orientation":"clockwise","parallel":0,"rotation":180.0,"section":"Tertiary section 0","winding":"Tertiary"}]} )");
        auto json_path_1009 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_inductance_bug_web_0_1009.json");
        std::ifstream json_file_1009(json_path_1009);
        json operatingPointData = json::parse(json_file_1009);

        Core core(coreData);
        OpenMagnetics::Coil winding(windingData);
        OperatingPoint operatingPoint(operatingPointData);
        MagnetizingInductance magnetizingInductance("ZHANG");
        double computedMagnetizingInductance = magnetizingInductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();
        // The computed magnetizing inductance must be finite and positive.
        CHECK(std::isfinite(computedMagnetizingInductance));
        CHECK(computedMagnetizingInductance > 0);
    }

    TEST_CASE("Test_Inductance_Bug_Web_1", "[physical-model][magnetizing-inductance][bug]") {
        settings.reset();
        clear_databases();

        json coreData = json::parse(R"({"distributorsInfo":null,"functionalDescription":{"coating":null,"gapping":[{"area":0.000315,"coordinates":[0,0,0],"distanceClosestNormalSurface":0.014498,"distanceClosestParallelSurface":0.011999999999999999,"length":0.000005,"sectionDimensions":[0.02,0.02],"shape":"round","type":"residual"},{"area":0.000164,"coordinates":[0.024563,0,0],"distanceClosestNormalSurface":0.014498,"distanceClosestParallelSurface":0.011999999999999999,"length":0.000005,"sectionDimensions":[0.005125,0.032],"shape":"irregular","type":"residual"},{"area":0.000164,"coordinates":[-0.024563,0,0],"distanceClosestNormalSurface":0.014498,"distanceClosestParallelSurface":0.011999999999999999,"length":0.000005,"sectionDimensions":[0.005125,0.032],"shape":"irregular","type":"residual"}],"material":"3C95","numberStacks":1,"shape":{"aliases":["EQ 50/32/20.0","EQ 50/20/32"],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.05},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.02},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.032},"D":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0145},"E":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.044},"F":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.02},"G":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.03205}},"family":"eq","familySubtype":null,"magneticCircuit":"open","name":"EQ 50/32/20","type":"standard"},"type":"twoPieceSet","magneticCircuit":"open"},"geometricalDescription":[{"coordinates":[0,0,0],"dimensions":null,"insulationMaterial":null,"machining":null,"material":"3C95","rotation":[3.141592653589793,3.141592653589793,0],"shape":{"aliases":["EQ 50/32/20.0","EQ 50/20/32"],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.05},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.02},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.032},"D":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0145},"E":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.044},"F":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.02},"G":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.03205}},"family":"eq","familySubtype":null,"magneticCircuit":"open","name":"EQ 50/32/20","type":"standard"},"type":"halfSet"},{"coordinates":[0,0,0],"dimensions":null,"insulationMaterial":null,"machining":null,"material":"3C95","rotation":[0,0,0],"shape":{"aliases":["EQ 50/32/20.0","EQ 50/20/32"],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.05},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.02},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.032},"D":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0145},"E":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.044},"F":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.02},"G":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.03205}},"family":"eq","familySubtype":null,"magneticCircuit":"open","name":"EQ 50/32/20","type":"standard"},"type":"halfSet"}],"manufacturerInfo":null,"name":"custom","processedDescription":{"columns":[{"area":0.000315,"coordinates":[0,0,0],"depth":0.02,"height":0.029,"minimumDepth":null,"minimumWidth":null,"shape":"round","type":"central","width":0.02},{"area":0.000164,"coordinates":[0.024563,0,0],"depth":0.032,"height":0.029,"minimumDepth":null,"minimumWidth":0.003001,"shape":"irregular","type":"lateral","width":0.005125},{"area":0.000164,"coordinates":[-0.024563,0,0],"depth":0.032,"height":0.029,"minimumDepth":null,"minimumWidth":0.003001,"shape":"irregular","type":"lateral","width":0.005125}],"depth":0.032,"effectiveParameters":{"effectiveArea":0.0003298035730425377,"effectiveLength":0.10383305467139949,"effectiveVolume":0.00003424451243054872,"minimumArea":0.0003141592653589793},"height":0.04,"thermalResistance":null,"width":0.05,"windingWindows":[{"angle":null,"area":0.00034799999999999995,"coordinates":[0.01,0],"height":0.029,"radialHeight":null,"sectionsAlignment":null,"sectionsOrientation":null,"shape":null,"width":0.011999999999999999}]}})");
        json windingData = json::parse(R"({"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.012082738036576609,"columnShape":"round","columnThickness":0.0020827380365766087,"columnWidth":0.012082738036576609,"coordinates":[0,0,0],"pins":null,"wallThickness":0.0017356150304805081,"windingWindows":[{"angle":null,"area":0.00025317549908941776,"coordinates":[0.017041369018288302,0,0],"height":0.025528769939038985,"radialHeight":null,"sectionsAlignment":null,"sectionsOrientation":null,"shape":"rectangular","width":0.00991726196342339}]}},"functionalDescription":[{"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":1,"wire":"Dummy"}],"layersDescription":null,"sectionsDescription":null,"turnsDescription":null})");
        json operatingPointData = json::parse(R"({"name":"Op. Point No. 1","conditions":{"ambientTemperature":100},"excitationsPerWinding":[{"frequency":50,"current":{"harmonics":{"amplitudes":[0,31.21831,0.348898,0.681642,11.946624,12.028013,10.143728,9.443625,10.985242,10.23761,10.161232,10.519392,11.473431,12.334118,15.790292,12.339356,3.782311,3.522273,3.435625,3.592056,3.338095,0.365192],"frequencies":[0,50,15000,25000,35000,45000,55000,65000,75000,85000,95000,105000,115000,125000,135000,145000,155000,165000,175000,185000,195000,205000]}},"voltage":{"harmonics":{"amplitudes":[0,1.136963,1.415776,2.520523,36.863994,43.911459,45.631752,50.314416,67.374707,71.317378,79.215135,90.682108,108.339425,126.882128,176.340424,143.759392,47.845538,47.502176,49.241268,54.317715,53.222847,6.05556],"frequencies":[0,50,15000,25000,35000,45000,55000,65000,75000,85000,95000,105000,115000,125000,135000,145000,155000,165000,175000,185000,195000,205000]}}}]})");

        Core core(coreData);
        OpenMagnetics::Coil winding(windingData);
        OperatingPoint operatingPoint(operatingPointData);
        MagnetizingInductance magnetizingInductance("ZHANG");
        double computedMagnetizingInductance = magnetizingInductance.calculate_inductance_from_number_turns_and_gapping(core, winding, &operatingPoint).get_magnetizing_inductance().get_nominal().value();
        // The computed magnetizing inductance must be finite and positive.
        CHECK(std::isfinite(computedMagnetizingInductance));
        CHECK(computedMagnetizingInductance > 0);
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

        auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "Error_inductance_with_Csv.json");
        auto mas = OpenMagneticsTesting::mas_loader(path.string());

        auto operatingPoint = mas.get_inputs().get_operating_points()[0];

        double magnetizingInductance =
            magnetizing_inductance.calculate_inductance_from_number_turns_and_gapping(mas.get_magnetic().get_core(), mas.get_magnetic().get_coil(), &operatingPoint).get_magnetizing_inductance().get_nominal().value();

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(magnetizingInductance, max_error * expectedValue));
    }


    // ABT #2074: a ground gap is one slab in the winding column, so the gap search may not step
    // past the column. It used to (1 mm steps, unbounded) and every step past it came back as a
    // core refusal, "gap of length 0.023005 m does not fit the winding column"; Henry saw ~5,000 of
    // them per corpus run. A target that no gap reaches is now refused once, with the inductance a
    // full-column gap still leaves. ETD 29 3C97, 666 turns: a full-column gap (22 mm) leaves 1.79 mH.
    TEST_CASE("Test_Gapping_Search_Refuses_A_Target_No_Gap_In_The_Column_Reaches", "[physical-model][magnetizing-inductance][gapping][abt-2074]") {
        settings.reset();
        clear_databases();
        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");
        prepare_test_parameters(0, 25, 20000, 666, 1e-3, {}, "ETD 29", "3C97", core, winding, inputs);

        try {
            magnetizing_inductance.calculate_gapping_from_number_turns_and_inductance(core, winding, &inputs, GappingType::GROUND);
            FAIL("a 1 mH target on 666 turns of an ETD 29 was given a gap");
        } catch (const GapException& e) {
            const std::string what = e.what();
            INFO(what);
            CHECK(what.find("No gap of core") != std::string::npos);
            CHECK(what.find("full") != std::string::npos);
            CHECK(what.find("does not fit the winding column") == std::string::npos);
        }
        settings.reset();
    }

    // ABT #2074: a target the column can just reach is solved inside the column, with the rounded gap
    // no longer than the column.
    TEST_CASE("Test_Gapping_Search_Near_The_Column_Limit_Stays_In_The_Column", "[physical-model][magnetizing-inductance][gapping][abt-2074]") {
        settings.reset();
        clear_databases();
        Core core;
        OpenMagnetics::Coil winding;
        OpenMagnetics::Inputs inputs;
        MagnetizingInductance magnetizing_inductance("ZHANG");
        prepare_test_parameters(0, 25, 20000, 666, 1.9e-3, {}, "ETD 29", "3C97", core, winding, inputs);
        const double columnHeight = core.get_processed_description()->get_columns()[0].get_height();

        auto gapping = magnetizing_inductance.calculate_gapping_from_number_turns_and_inductance(core, winding, &inputs, GappingType::GROUND);
        INFO("gap " << gapping[0].get_length() * 1e3 << " mm in a " << columnHeight * 1e3 << " mm column");
        CHECK(gapping[0].get_length() > 15e-3);
        CHECK(gapping[0].get_length() <= columnHeight);
        settings.reset();
    }

}  // namespace

// ABT #331: open-core (drum) magnetizing inductance, validated against PUBLISHED vendor data
// rather than another MKF result. Fair-Rite sells bare drum cores ("bobbins") with AL printed on
// every part page; the four parts below are the ones whose dimension-letter mapping was
// WEIGHT-VERIFIED (computed volume x material density reproduces the listed grams to 3-4%).
//
// The model is the demagnetising-factor bracket (see MagnetizingInductance.cpp): upper bound =
// flange-envelope spheroid, lower bound = envelope air return in series with the ferrite post,
// estimate = log-midpoint. Measured envelope on this exact set: mean 5.4%, max 9.9%. The 12%
// tolerance below is that envelope plus margin; if this test starts failing the MODEL drifted,
// so investigate rather than widen.
//
// Note the physics this pins: Fair-Rite publishes the SAME AL for 43-material (mu_i 800) and
// 77-material (mu_i 2000) variants of one geometry, because an open core's inductance saturates
// at the geometry-set limit ~1/N_d. The model must reproduce that material-insensitivity.
TEST_CASE("Test_Open_Core_Drum_Inductance_Matches_FairRite_AL", "[physical-model][magnetizing-inductance][drum][open-core]") {
    settings.reset();
    clear_databases();

    struct Reference {
        std::string shapeName;
        std::string materialName;
        double testCoilTurns;
        double publishedAlNanoHenry;   // Fair-Rite part page, 1 kHz < 10 gauss
    };
    std::vector<Reference> references = {
        {"Bobbin 9643001015", "43", 75, 38.0},
        {"Bobbin 9677282509", "77", 55, 95.0},
        {"Bobbin 9677182209", "77", 95, 65.0},
        {"Bobbin 9677282009", "77", 40, 100.0},
    };

    for (const auto& reference : references) {
        auto core = OpenMagneticsTesting::get_quick_core(reference.shapeName, json::array(), 1, reference.materialName);
        REQUIRE(core.get_functional_description().get_type() == CoreType::OPEN_SHAPE);

        double inductance = MagnetizingInductance::calculate_open_core_magnetizing_inductance(
            core, reference.testCoilTurns, 25);
        double publishedInductance = reference.publishedAlNanoHenry * 1e-9 * pow(reference.testCoilTurns, 2);
        UNSCOPED_INFO(reference.shapeName << ": model " << inductance * 1e6 << " uH vs published "
                      << publishedInductance * 1e6 << " uH ("
                      << (inductance - publishedInductance) / publishedInductance * 100 << "%)");
        CHECK_THAT(inductance, Catch::Matchers::WithinRel(publishedInductance, 0.12));
    }

    // Gapping an open core is meaningless and must throw, not silently compute.
    auto gapped = json::array();
    gapped.push_back(json{{"type", "additive"}, {"length", 0.0005}});
    auto gappedCore = OpenMagneticsTesting::get_quick_core("Bobbin 9643001015", gapped, 1, "43");
    CHECK_THROWS(MagnetizingInductance::calculate_open_core_magnetizing_inductance(gappedCore, 10, 25));
    settings.reset();
}

// ABT #331: end-to-end DEMO on a real commercial unshielded drum inductor — WE-TI 7447720470
// (Wurth, 47 uH +-10%, DCR typ 89 mOhm, envelope flanges 7.8/5.0 mm). The vendor does not publish
// turns or internal core dims, so the fixture is a RECONSTRUCTION, documented as such: the
// open-core model's AL is post-diameter-insensitive (envelope-dominated), giving N = 48; the wire
// that then reproduces the vendor DCR (0.42 mm Cu, 4 layers) also FITS the groove. Three vendor
// numbers, two assumptions, closed loop. This is a consistency demo on a finished part —
// the strict model validation is the Fair-Rite AL test above.
TEST_CASE("Test_Open_Core_WE_TI_Reconstruction_Consistency", "[physical-model][magnetizing-inductance][drum][open-core]") {
    settings.reset();
    clear_databases();
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "we_ti_7447720470_reconstructed.json");
    std::ifstream file(path);
    REQUIRE(file.good());
    auto masJson = nlohmann::json::parse(file);
    OpenMagnetics::Core core(masJson["magnetic"]["core"]);
    core.process_data();
    REQUIRE(core.get_functional_description().get_type() == CoreType::OPEN_SHAPE);

    double numberTurns = masJson["magnetic"]["coil"]["functionalDescription"][0]["numberTurns"];
    double inductance = MagnetizingInductance::calculate_open_core_magnetizing_inductance(core, numberTurns, 25);
    // Vendor L = 47 uH +-10%; the model must land inside vendor tolerance + its own 12% envelope.
    UNSCOPED_INFO("model L = " << inductance * 1e6 << " uH vs vendor 47 uH +-10%");
    CHECK_THAT(inductance, Catch::Matchers::WithinRel(47e-6, 0.20));

    // DCR of the reconstructed winding: N x MLT(post + buildup) x rho/A must reproduce the
    // vendor's typ 89 mOhm within winding-tolerance slack.
    double wireDiameter = 0.00042;
    double wireArea = std::numbers::pi / 4 * pow(wireDiameter, 2);
    double buildup = 4 * 0.000458;   // 4 layers of 0.458 OD
    double meanTurnDiameter = 0.0035 + buildup;
    double length = numberTurns * std::numbers::pi * meanTurnDiameter;
    double dcr = 1.72e-8 * length / wireArea;
    UNSCOPED_INFO("reconstructed DCR = " << dcr * 1000 << " mOhm vs vendor typ 89");
    CHECK_THAT(dcr, Catch::Matchers::WithinRel(0.089, 0.30));
    settings.reset();
}

// ABT #366: shielded drum (drumRing). No vendor publishes AL for an assembled drum+ring pair
// (verified industry-wide 2026-07-29: Ferroxcube/Fair-Rite publish bare-drum AL only, with a
// defined winding; ACME publishes matched DR+SRI geometry but no electricals), so this pins the
// PHYSICS BRACKETS instead: the closed-through-clearance-gaps circuit must yield strictly MORE
// inductance than the bare drum (the ferrite ring replaces most of the air return) and strictly
// LESS than the same ferrite circuit with the clearance gaps shorted (core-only reluctance).
// The analytic two-annular-gap estimate anchors the magnitude: for a mu_i ~2000 ferrite the two
// clearance gaps dominate the closed circuit, and fringing can only lower their reluctance.
// FEM validation of the absolute value is the follow-up recorded in ABT #366.
TEST_CASE("Test_Drum_Ring_Inductance_Brackets", "[physical-model][magnetizing-inductance][drum-ring]") {
    settings.reset();
    clear_databases();
    double numberTurns = 20;

    // Bare drum with the SAME drum dimensions (custom shape: no bare-drum MAS record exists
    // for the ACME DR2.3 drum, only the paired drumRing record).
    json bareShapeJson = {
        {"magneticCircuit", "open"}, {"type", "custom"}, {"family", "drum"},
        {"aliases", json::array()}, {"name", "DR 2.3 bare"},
        {"dimensions", {
            {"A", {{"nominal", 0.0023}}}, {"B", {{"nominal", 0.001}}}, {"C", {{"nominal", 0.0011}}},
            {"D", {{"nominal", 0.00021}}}, {"E", {{"nominal", 0.00058}}}, {"F", {{"nominal", 0.00021}}}}}
    };
    json bareCoreJson;
    bareCoreJson["functionalDescription"] = {
        {"type", "openShape"}, {"material", "3C90"}, {"shape", bareShapeJson},
        {"gapping", json::array()}, {"numberStacks", 1}};
    Core bareCore(bareCoreJson);
    bareCore.process_data();
    double bareInductance = MagnetizingInductance::calculate_open_core_magnetizing_inductance(bareCore, numberTurns, 25);

    // The shielded assembly, from the MAS drumRing record (same drum + SRI 3x2.4x1.05 ring).
    auto core = OpenMagneticsTesting::get_quick_core("DR 2.3 + SRI 3.0", json::array(), 1, "3C90");
    json windingData = json::parse(
        R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary",
            "numberParallels": 1, "numberTurns": 20, "wire": "Dummy"}]})");
    OpenMagnetics::Coil winding(windingData);
    MagnetizingInductance magnetizingInductanceModel("ZHANG");
    double ringInductance = magnetizingInductanceModel
        .calculate_inductance_from_number_turns_and_gapping(core, winding, nullptr)
        .get_magnetizing_inductance().get_nominal().value();

    // Upper limit: the same ferrite circuit with the clearance gaps shorted
    // (get_core_reluctance() is core PLUS gaps; the ungapped field is the core alone).
    auto reluctanceModel = ReluctanceModel::factory(ReluctanceModels::ZHANG);
    auto reluctanceOutput = reluctanceModel->get_core_reluctance(core);
    double coreOnlyReluctance = reluctanceOutput.get_ungapped_core_reluctance().value();
    double shortedGapsInductance = pow(numberTurns, 2) / coreOnlyReluctance;

    UNSCOPED_INFO("bare " << bareInductance * 1e6 << " uH < ring " << ringInductance * 1e6
                          << " uH < shorted-gaps " << shortedGapsInductance * 1e6 << " uH");
    CHECK(std::isfinite(ringInductance));
    CHECK(bareInductance < ringInductance);
    CHECK(ringInductance < shortedGapsInductance);

    // The two annular gaps must be evaluated INDIVIDUALLY, combine in SERIES (single column:
    // both land on the wound post, never in the parallel lateral term), and ENGAGE the fringing
    // machinery (factor > 1 lowers their reluctance, raising L over the sharp-gap value).
    REQUIRE(reluctanceOutput.get_reluctance_per_gap().has_value());
    auto reluctancePerGap = reluctanceOutput.get_reluctance_per_gap().value();
    REQUIRE(reluctancePerGap.size() == 2);
    double gappingReluctance = reluctanceOutput.get_gapping_reluctance().value();
    CHECK_THAT(gappingReluctance,
               Catch::Matchers::WithinRel(reluctancePerGap[0].get_reluctance() + reluctancePerGap[1].get_reluctance(), 1e-9));
    CHECK(reluctanceOutput.get_maximum_fringing_factor().value() > 1.0);

    // Analytic anchor: two annular clearance gaps in series, no fringing.
    double gapLength = (0.0024 - 0.0023) / 2;
    double gapArea = 2 * std::numbers::pi * ((0.0023 + 0.0024) / 4) * 0.00021;
    double basicTwoGapReluctance = 2 * gapLength / (4e-7 * std::numbers::pi * gapArea);
    double basicEstimate = pow(numberTurns, 2) / basicTwoGapReluctance;
    UNSCOPED_INFO("analytic two-gap estimate " << basicEstimate * 1e6 << " uH");
    CHECK(ringInductance > 0.5 * basicEstimate);
    CHECK(ringInductance < 3.0 * basicEstimate);
    settings.reset();
}

// ABT #417: drumRing's two structural annular-clearance gaps are DERIVED (Core::process_gap
// synthesizes them from A/K/D/F) — functionalDescription.gapping is legitimately empty in the
// raw MAS record, unlike a normal gapped E-core where the user/adviser always specifies real
// entries. Core's free from_json (used whenever a Core is deserialized as a struct MEMBER —
// e.g. Magnetic::from_json, which every PyOM/WASM binding reaches for a loaded MAS document)
// never calls process_data()/process_gap(), unlike the Core(json) constructor. A drumRing core
// arriving via that path kept gapping==[], and ReluctanceModel::get_gapping_reluctance treated
// "empty" as "no gap" instead of "not yet derived" — silently dropping the dominant reluctance
// term and reporting an inductance 3.8-10.7x too high, so calculate_saturation_current
// under-reported Isat by the same factor. Fixed in calculate_inductance_and_magnetic_flux_density:
// self-heals (calls core.process_gap()) whenever a drumRing core arrives with empty gapping.
TEST_CASE("Test_ABT417_DrumRing_Saturation_Current_Matches_Model_Inductance_Unprocessed", "[physical-model][magnetizing-inductance][drum-ring][bug]") {
    settings.reset();
    clear_databases();
    int64_t numberTurns = 20;

    // Raw catalog-style json (shape as a bare name string, gapping empty) round-tripped through
    // Core's free from_json — reproduces the state a Magnetic loaded via Magnetic::from_json
    // arrives in (NOT the Core(json) constructor, which self-processes and would mask this bug).
    json rawCoreJson;
    rawCoreJson["functionalDescription"]["name"] = "ABT417Test";
    rawCoreJson["functionalDescription"]["type"] = "pieceAndPlate";
    rawCoreJson["functionalDescription"]["material"] = "3C90";
    rawCoreJson["functionalDescription"]["shape"] = "DR 2.3 + SRI 3.0";
    rawCoreJson["functionalDescription"]["gapping"] = json::array();
    rawCoreJson["functionalDescription"]["numberStacks"] = 1;

    Core unprocessedCore;
    OpenMagnetics::from_json(rawCoreJson, unprocessedCore);
    REQUIRE(unprocessedCore.get_functional_description().get_gapping().empty());
    // Populate processedDescription (effective area etc.) WITHOUT deriving the structural
    // clearance gaps — process_data() and process_gap() are separate calls, and this is the
    // partially-processed state a core can plausibly reach outside the Core(json) constructor
    // (which always runs both together). get_effective_area() requires processedDescription,
    // so calculate_saturation_current would throw CoreNotProcessedException without this —
    // masking the numeric bug behind an exception instead of reproducing the ticket's reported
    // finite-but-wrong Isat.
    unprocessedCore.process_data();
    REQUIRE(unprocessedCore.get_functional_description().get_gapping().empty());

    json windingData = json::parse(
        R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary",
            "numberParallels": 1, "numberTurns": 20, "wire": "Dummy"}]})");
    OpenMagnetics::Coil winding(windingData);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(unprocessedCore);
    magnetic.set_coil(winding);

    MagnetizingInductance magnetizingInductanceModel("ZHANG");
    double modelInductance = magnetizingInductanceModel
        .calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil(), nullptr)
        .get_magnetizing_inductance().get_nominal().value();

    double saturationCurrent = magnetic.calculate_saturation_current(25, false);
    double bSat = magnetic.get_mutable_core().get_magnetic_flux_density_saturation(25, false);
    double effectiveArea = magnetic.get_mutable_core().get_effective_area();
    double impliedInductance = bSat * numberTurns * effectiveArea / saturationCurrent;

    // Ground truth: the SAME shape/material/turns, but the core built through the
    // Core(json) constructor path (which always derives the clearance gaps regardless of
    // this fix — see Test_Drum_Ring_Inductance_Brackets). Comparing modelInductance only
    // against impliedInductance is not enough: before the fix, BOTH go through the same
    // unhealed generic path and agree with each other while being 30x too high together.
    // This reference catches that — it must independently derive the same, correct answer.
    auto referenceCore = OpenMagneticsTesting::get_quick_core("DR 2.3 + SRI 3.0", json::array(), 1, "3C90");
    double referenceInductance = magnetizingInductanceModel
        .calculate_inductance_from_number_turns_and_gapping(referenceCore, winding, nullptr)
        .get_magnetizing_inductance().get_nominal().value();

    UNSCOPED_INFO("model inductance " << modelInductance * 1e6 << " uH, saturation-implied inductance "
                                      << impliedInductance * 1e6 << " uH, reference (always-processed) "
                                      << referenceInductance * 1e6 << " uH");
    CHECK_THAT(impliedInductance, Catch::Matchers::WithinRel(modelInductance, 0.05));
    CHECK_THAT(modelInductance, Catch::Matchers::WithinRel(referenceInductance, 0.05));

    CHECK(unprocessedCore.get_shape_family() == CoreShapeFamily::DRUM_RING);
    settings.reset();
}

// ABT #357: molded composite body (WE-MAPI class). The whole circuit is IN the low-mu
// material — no discrete gap dilutes it — so L must scale essentially linearly with the
// composite permeability (Kool Mu 60 vs 26 => ratio ~2.31), and the magnitude must agree
// with mu0 * mu * N^2 * Ae / le computed from the piece's own effective parameters. This is
// the property the ABT #357 phase-2 material extraction inverts (mu from measured L0).
TEST_CASE("Test_Molded_Inductance_Permeability_Scaling", "[physical-model][magnetizing-inductance][molded]") {
    settings.reset();
    clear_databases();
    double numberTurns = 10;
    json windingData = json::parse(
        R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary",
            "numberParallels": 1, "numberTurns": 10, "wire": "Dummy"}]})");

    auto buildCore = [](const std::string& materialName) {
        json shapeJson = {
            {"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "molded"},
            {"aliases", json::array()}, {"name", "MAPI-like 4020"},
            {"dimensions", {
                {"A", {{"nominal", 0.0041}}}, {"B", {{"nominal", 0.0021}}}, {"C", {{"nominal", 0.0041}}},
                {"D", {{"nominal", 0.0014}}}, {"E", {{"nominal", 0.0030}}}, {"F", {{"nominal", 0.0012}}}}}
        };
        json coreJson;
        coreJson["functionalDescription"] = {
            {"type", "closedShape"}, {"material", materialName}, {"shape", shapeJson},
            {"gapping", json::array()}, {"numberStacks", 1}};
        Core core(coreJson);
        core.process_data();
        core.process_gap();
        return core;
    };

    MagnetizingInductance magnetizingInductanceModel("ZHANG");
    auto core26 = buildCore("Kool Mµ 26");
    auto core60 = buildCore("Kool Mµ 60");
    double inductance26 = magnetizingInductanceModel
        .calculate_inductance_from_number_turns_and_gapping(core26, OpenMagnetics::Coil(windingData))
        .get_magnetizing_inductance().get_nominal().value();
    double inductance60 = magnetizingInductanceModel
        .calculate_inductance_from_number_turns_and_gapping(core60, OpenMagnetics::Coil(windingData))
        .get_magnetizing_inductance().get_nominal().value();

    UNSCOPED_INFO("L(mu26) = " << inductance26 * 1e9 << " nH, L(mu60) = " << inductance60 * 1e9 << " nH");
    CHECK(std::isfinite(inductance26));
    CHECK(inductance26 > 0);
    // Linear-in-mu within the tolerance the vendor DC-bias fits leave at H=0.
    CHECK(inductance60 / inductance26 > 1.9);
    CHECK(inductance60 / inductance26 < 2.5);

    // Magnitude agrees with the closed-circuit formula on the piece's own effective
    // parameters (same le/Ae source, so the band only absorbs the permeability fit at H=0).
    auto effectiveParameters = core26.get_processed_description().value().get_effective_parameters();
    double handEstimate = 4e-7 * std::numbers::pi * 26 * pow(numberTurns, 2) *
                          effectiveParameters.get_effective_area() / effectiveParameters.get_effective_length();
    CHECK(inductance26 > 0.7 * handEstimate);
    CHECK(inductance26 < 1.4 * handEstimate);
    settings.reset();
}

// ABT #357 phase 2: MKF's molded model against the REAL WE-MAPI catalogue. The reconstruction
// (scripts/fit_we_mapi_internals.py) fitted each part's cavity geometry, turn count, wire gauge
// and composite mu_eff0 from PUBLIC data only — datasheet L0 + DCR + body dimensions + the
// measured L(I) curves — using a PYTHON MIRROR of CorePieceMolded's c1/c2. That mirror is the
// risk this test exists to kill: if the C++ and the fit ever disagree, every reconstruction
// silently rots. So here MKF itself rebuilds each part from its reconstructed letters and must
// reproduce the part's DATASHEET inductance.
//
// Honest scope (see the fit's own degeneracy analysis, carried in we_mapi_fitted_materials.json):
// the public data pins N, wire gauge and the mu(H) rolloff shape well, but mu_eff0 only to about
// +-25% (L0/DCR/packing leave a broad valley; a B-at-drop physical band selects the branch). So
// the assertion is on MKF-vs-DATASHEET agreement at the fit's own residual level, NOT on mu being
// the true composite permeability.
TEST_CASE("Test_Molded_WE_MAPI_Reconstructions_Reproduce_Datasheet_Inductance",
          "[physical-model][magnetizing-inductance][molded][we-mapi]") {
    settings.reset();
    clear_databases();
    namespace fs = std::filesystem;
    auto reconstructionsPath = fs::path{std::source_location::current().file_name()}
                                   .parent_path().append("testData").append("we_mapi_reconstructions.json");
    auto stubsPath = fs::path{std::source_location::current().file_name()}
                         .parent_path().append("testData").append("we_mapi_datasheet_stubs.ndjson");
    std::ifstream reconstructionsFile(reconstructionsPath);
    std::ifstream stubsFile(stubsPath);
    REQUIRE(reconstructionsFile.good());
    REQUIRE(stubsFile.good());
    json reconstructions = json::parse(reconstructionsFile);
    REQUIRE(reconstructions.size() == 183);

    // Datasheet inductance per order code, straight from the RedExpert-derived stubs.
    std::map<std::string, double> datasheetInductance;
    std::string stubLine;
    while (std::getline(stubsFile, stubLine)) {
        if (stubLine.empty()) continue;
        json stub = json::parse(stubLine);
        auto& manufacturerInfo = stub.at("manufacturerInfo");
        datasheetInductance[manufacturerInfo.at("reference").get<std::string>()] =
            manufacturerInfo.at("datasheetInfo").at("electrical").at(0)
                            .at("inductance").at("nominal").get<double>();
    }
    REQUIRE(datasheetInductance.size() == 183);

    MagnetizingInductance magnetizingInductanceModel("ZHANG");
    std::vector<double> relativeErrors;
    size_t evaluated = 0;
    size_t withinTwentyPercent = 0;
    for (auto& reconstruction : reconstructions) {
        auto orderCode = reconstruction.at("orderCode").get<std::string>();
        // The 4012 family is a 1-2 turn (almost certainly stamped/flat-wire) construction the
        // round-wire reconstruction cannot represent — the fit flags it, so it is excluded here
        // rather than silently widening the band for everyone else.
        if (reconstruction.at("caseCode").get<std::string>() == "4012") continue;
        double numberTurns = reconstruction.at("numberTurns").get<double>();
        double permeability = reconstruction.at("mu_eff0").get<double>();

        // Inline material carrying the FITTED composite permeability (labelled as such; this is
        // a fixture, not a catalogue material record).
        json materialJson = {
            {"name", "WE metal alloy (fitted mu_eff0)"}, {"type", "custom"},
            {"material", "powder"}, {"materialComposition", "proprietary"},
            {"manufacturerInfo", {{"name", "FITTED — not vendor material data"}}},
            {"permeability", {{"initial", {{"value", permeability}}}}},
            {"resistivity", json::array({{{"value", 1.0}}})},
            {"density", 5500},
            {"curieTemperature", 500},
            {"saturation", json::array({{{"magneticFluxDensity", 1.2}, {"magneticField", 40000}, {"temperature", 25}}})},
            {"volumetricLosses", {{"default", json::array({{{"method", "lossFactor"},
                {"factors", json::array({{{"value", 1e-5}, {"frequency", 100000}}})}}})}}}
        };
        json shapeJson = {
            {"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "molded"},
            {"aliases", json::array()}, {"name", "WE-MAPI " + orderCode},
            {"dimensions", {
                {"A", {{"nominal", reconstruction.at("A").get<double>()}}},
                {"B", {{"nominal", reconstruction.at("B").get<double>()}}},
                {"C", {{"nominal", reconstruction.at("C").get<double>()}}},
                {"D", {{"nominal", reconstruction.at("D").get<double>()}}},
                {"E", {{"nominal", reconstruction.at("E").get<double>()}}},
                {"F", {{"nominal", reconstruction.at("F").get<double>()}}}}}
        };
        json coreJson;
        coreJson["functionalDescription"] = {
            {"type", "closedShape"}, {"material", materialJson}, {"shape", shapeJson},
            {"gapping", json::array()}, {"numberStacks", 1}};
        Core core(coreJson);
        core.process_data();
        core.process_gap();

        json coilJson;
        coilJson["bobbin"] = "Dummy";
        coilJson["functionalDescription"] = json::array({{
            {"name", "winding 0"}, {"numberTurns", int(numberTurns)}, {"numberParallels", 1},
            {"isolationSide", "primary"}, {"wire", "Dummy"}}});
        OpenMagnetics::Coil coil(coilJson);

        double modelInductance = magnetizingInductanceModel
            .calculate_inductance_from_number_turns_and_gapping(core, coil)
            .get_magnetizing_inductance().get_nominal().value();
        double relativeError = modelInductance / datasheetInductance.at(orderCode) - 1;
        relativeErrors.push_back(std::abs(relativeError));
        if (std::abs(relativeError) < 0.20) withinTwentyPercent++;
        evaluated++;
    }
    REQUIRE(evaluated >= 175);

    std::sort(relativeErrors.begin(), relativeErrors.end());
    double medianRelativeError = relativeErrors[relativeErrors.size() / 2];
    double worstRelativeError = relativeErrors.back();
    UNSCOPED_INFO("MKF vs datasheet over " << evaluated << " WE-MAPI parts: median |err| "
                  << medianRelativeError * 100 << "%, worst " << worstRelativeError * 100
                  << "%, within 20%: " << withinTwentyPercent);
    // The fit's own median |L residual| is 3.9%; MKF must land in the same neighbourhood, which
    // is what proves the C++ and the python mirror agree. A drift here means CorePieceMolded's
    // sectioning changed and the reconstructions need refitting — do NOT widen this band.
    CHECK(medianRelativeError < 0.10);
    CHECK(withinTwentyPercent > evaluated * 8 / 10);
    settings.reset();
}

// ABT #362: semi-shielded drum — mixed-material sectioned reluctance (ferrite drum + magnetic-
// epoxy shell). No public vendor data pairs a drum+glue geometry with L (the 290-part
// validation set is confidential and runs on heimdall's side), so this pins the PHYSICS:
// the glue return must beat the bare drum's air return; a higher glue mu must give more L
// (monotonicity, saturating towards the all-ferrite limit); and every missing material link
// must THROW (no-fallbacks: the model never assumes air or guesses a mu). Glue stand-ins are
// existing catalogue powder/ferrite materials — real glue records are ABT #364.
TEST_CASE("Test_Drum_Semishielded_Inductance", "[physical-model][magnetizing-inductance][drum-semishielded]") {
    settings.reset();
    clear_databases();
    double numberTurns = 20;

    auto buildCore = [](json coatingJson) {
        json shapeJson = {
            {"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "drumSemishielded"},
            {"aliases", json::array()}, {"name", "LQS-like 4018"},
            {"dimensions", {
                {"A", {{"nominal", 0.0038}}}, {"B", {{"nominal", 0.0018}}}, {"C", {{"nominal", 0.0015}}},
                {"D", {{"nominal", 0.0004}}}, {"E", {{"nominal", 0.0010}}}, {"F", {{"nominal", 0.0004}}},
                {"J", {{"nominal", 0.0040}}}, {"K", {{"nominal", 0.0040}}}, {"L", {{"nominal", 0.0018}}}}}
        };
        json coreJson;
        coreJson["functionalDescription"] = {
            {"type", "pieceAndPlate"}, {"material", "3C90"}, {"shape", shapeJson},
            {"gapping", json::array()}, {"numberStacks", 1}};
        if (!coatingJson.is_null()) {
            coreJson["functionalDescription"]["coating"] = coatingJson;
        }
        Core core(coreJson);
        core.process_data();
        core.process_gap();
        return core;
    };
    auto glueCoating = [](const std::string& materialName) {
        return json{{"type", "magneticEpoxy"}, {"thickness", {{"nominal", 0.0001}}}, {"material", materialName}};
    };

    // Bare drum with the SAME drum dimensions: the air return the glue replaces.
    json bareShapeJson = {
        {"magneticCircuit", "open"}, {"type", "custom"}, {"family", "drum"},
        {"aliases", json::array()}, {"name", "LQS drum bare"},
        {"dimensions", {
            {"A", {{"nominal", 0.0038}}}, {"B", {{"nominal", 0.0018}}}, {"C", {{"nominal", 0.0015}}},
            {"D", {{"nominal", 0.0004}}}, {"E", {{"nominal", 0.0010}}}, {"F", {{"nominal", 0.0004}}}}}
    };
    json bareCoreJson;
    bareCoreJson["functionalDescription"] = {
        {"type", "openShape"}, {"material", "3C90"}, {"shape", bareShapeJson},
        {"gapping", json::array()}, {"numberStacks", 1}};
    Core bareCore(bareCoreJson);
    bareCore.process_data();
    double bareInductance = MagnetizingInductance::calculate_open_core_magnetizing_inductance(bareCore, numberTurns, 25);

    double glue26Inductance = MagnetizingInductance::calculate_semishielded_drum_magnetizing_inductance(
        buildCore(glueCoating("Kool Mµ 26")), numberTurns, 25);
    double glue60Inductance = MagnetizingInductance::calculate_semishielded_drum_magnetizing_inductance(
        buildCore(glueCoating("Kool Mµ 60")), numberTurns, 25);
    // All-ferrite shell = the drumRing-like upper limit of the same circuit.
    double ferriteShellInductance = MagnetizingInductance::calculate_semishielded_drum_magnetizing_inductance(
        buildCore(glueCoating("3C90")), numberTurns, 25);

    UNSCOPED_INFO("bare " << bareInductance * 1e6 << " uH < glue26 " << glue26Inductance * 1e6
                          << " uH < glue60 " << glue60Inductance * 1e6 << " uH < ferrite shell "
                          << ferriteShellInductance * 1e6 << " uH");
    CHECK(std::isfinite(glue26Inductance));
    CHECK(bareInductance < glue26Inductance);
    CHECK(glue26Inductance < glue60Inductance);
    CHECK(glue60Inductance < ferriteShellInductance);

    // The routed path (family gate) must agree with the direct call.
    json windingData = json::parse(
        R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary",
            "numberParallels": 1, "numberTurns": 20, "wire": "Dummy"}]})");
    MagnetizingInductance magnetizingInductanceModel("ZHANG");
    double routedInductance = magnetizingInductanceModel
        .calculate_inductance_from_number_turns_and_gapping(buildCore(glueCoating("Kool Mµ 26")),
                                                            OpenMagnetics::Coil(windingData))
        .get_magnetizing_inductance().get_nominal().value();
    CHECK_THAT(routedInductance, Catch::Matchers::WithinRel(glue26Inductance, 1e-9));

    // No-fallbacks: every missing link throws.
    CHECK_THROWS(MagnetizingInductance::calculate_semishielded_drum_magnetizing_inductance(
        buildCore(json()), numberTurns, 25));  // no coating at all
    CHECK_THROWS(MagnetizingInductance::calculate_semishielded_drum_magnetizing_inductance(
        buildCore(json{{"type", "epoxy"}, {"thickness", {{"nominal", 0.0001}}}}), numberTurns, 25));  // non-magnetic coating
    CHECK_THROWS(MagnetizingInductance::calculate_semishielded_drum_magnetizing_inductance(
        buildCore(glueCoating("No Such Glue")), numberTurns, 25));  // unknown shell material
    settings.reset();
}

TEST_CASE("Test_ABT635_Gapping_From_Inductance_Without_ProcessedDescription",
          "[physical-model][magnetizing-inductance][bug]") {
    // ABT #635: calculate_gapping_from_number_turns_and_inductance() dereferenced
    // core.get_processed_description() unconditionally, so a Core built from a bare (and perfectly
    // legal) functionalDescription -- no processedDescription -- failed with an opaque
    // "bad optional access" instead of being processed. Every neighbouring entry point accepts
    // that same core, so the API was inconsistent as well as unhelpful.
    settings.reset();
    clear_databases();

    // identical core TWICE: once with only a functionalDescription, once pre-processed.
    const std::string coreJson = R"({"name": "abt635", "functionalDescription": {
        "type": "twoPieceSet", "material": "3C95",
        "shape": {"type": "custom", "family": "pq", "name": "abt635 pq",
                  "dimensions": {"A": 0.0273, "B": 0.009, "C": 0.019, "D": 0.0057,
                                 "E": 0.0225, "F": 0.012, "G": 0.0155}},
        "gapping": [], "numberStacks": 1}})";

    Core bare = json::parse(coreJson);
    Core processed = json::parse(coreJson);
    processed.process_data();
    processed.process_gap();

    OpenMagnetics::Coil coil = json::parse(R"({"bobbin": "Dummy", "functionalDescription": [
        {"name": "PRI", "numberTurns": 10, "numberParallels": 1, "isolationSide": "primary",
         "wire": "Round 21.0 - Heavy Build"}]})");
    OpenMagnetics::Inputs inputs = json::parse(R"({"designRequirements": {
        "name": "abt635", "magnetizingInductance": {"nominal": 3.3e-06}, "turnsRatios": []},
        "operatingPoints": []})");

    MagnetizingInductance model(Defaults().reluctanceModelDefault);
    GappingType gappingType = magic_enum::enum_cast<GappingType>("GROUND").value();

    // the ticket's exact repro: this used to throw
    std::vector<CoreGap> fromBare;
    REQUIRE_NOTHROW(fromBare = model.calculate_gapping_from_number_turns_and_inductance(
        bare, coil, &inputs, gappingType, 6));
    REQUIRE(fromBare.size() > 0);
    REQUIRE(fromBare[0].get_length() > 0);

    // and it must agree with the pre-processed core -- processing internally must not change the
    // answer, only remove the need for the caller to have done it.
    auto fromProcessed = model.calculate_gapping_from_number_turns_and_inductance(
        processed, coil, &inputs, gappingType, 6);
    REQUIRE(fromProcessed.size() == fromBare.size());
    CHECK_THAT(fromBare[0].get_length(),
               Catch::Matchers::WithinRel(fromProcessed[0].get_length(), 1e-9));

    settings.reset();
}

// ABT #907: on a CCM flybuck the magnetizing-current DC anchor was taken as
// peak(I1) - pkpk(I1)/2, but the primary's pk-pk excursion contains the
// commutation step (the secondary's ampere-turns handed over at turn-on),
// which moves no flux. That understated the DC anchor (0.420 A instead of
// 0.583 A) and with it B_dc (0.191 T instead of 0.266 T) and B_peak (0.259 T
// instead of 0.334 T) — the saturation margin read 28% better than reality.
// The anchor is now peak-based: peak(i_mag) == peak(I1) (at the primary peak
// the secondary is off, so the primary carries the whole magnetizing MMF),
// minus the magnetizing excursion above its mean from the volt-second
// integral. User-reported design: High Flux 160 T10/6/5, 56:28 turns,
// 170 kHz flybuck. The AC flux swing (67.8 mT peak) is pinned by Faraday's
// law on the primary volt-seconds and must not move with this fix.
TEST_CASE("Test_ABT907_Flybuck_Flux_Density_Dc_Offset_Uses_Net_Mmf_Anchor",
          "[physical-model][magnetizing-inductance][flyback][abt907]") {
    settings.reset();
    clear_databases();
    namespace fs = std::filesystem;
    auto masPath = fs::path{std::source_location::current().file_name()}
                       .parent_path().append("testData").append("flybuck_hf160_250u_abt907.json");
    std::ifstream masFile(masPath);
    REQUIRE(masFile.good());
    json masJson = json::parse(masFile);

    OpenMagnetics::Magnetic magnetic(masJson["magnetic"]);
    OperatingPoint operatingPoint(masJson["inputs"]["operatingPoints"][0]);

    MagnetizingInductance magnetizingInductanceModel;
    auto [inductanceOutput, fluxDensity] = magnetizingInductanceModel.calculate_inductance_and_magnetic_flux_density(
        magnetic.get_core(), magnetic.get_coil(), &operatingPoint);
    auto fluxProcessed = fluxDensity.get_processed().value();

    // AC part: Faraday on the primary voltage waveform gives 75.9 uVs over
    // 56 turns x 10 mm2 -> 135.6 mT pk-pk. Model-independent; unchanged by the fix.
    CHECK_THAT(fluxProcessed.get_peak_to_peak().value(), Catch::Matchers::WithinRel(0.1356, 0.02));

    // DC part: the magnetizing current i_m = I1 + (N2/N1)*I2 is a clean triangle
    // running 0.4347 -> 0.7326 A (verified against the pointwise volt-second
    // integral); its midpoint 0.5837 A through the core reluctance gives 0.266 T.
    CHECK_THAT(fluxProcessed.get_offset(), Catch::Matchers::WithinRel(0.2657, 0.02));
    CHECK_THAT(fluxProcessed.get_peak().value(), Catch::Matchers::WithinRel(0.3335, 0.02));

    // The anchor itself: peak(i_mag) must coincide with the measured primary peak
    // (0.7326 A), and the valley with the commutation-corrected 0.4347 A — NOT the
    // primary waveform's own 0.108 A minimum.
    auto magnetizingCurrentWaveform =
        operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current().value().get_waveform().value();
    const std::vector<double>& magnetizingCurrentData = magnetizingCurrentWaveform.get_data();
    double magnetizingCurrentMaximum = *max_element(magnetizingCurrentData.begin(), magnetizingCurrentData.end());
    double magnetizingCurrentMinimum = *min_element(magnetizingCurrentData.begin(), magnetizingCurrentData.end());
    CHECK_THAT(magnetizingCurrentMaximum, Catch::Matchers::WithinRel(0.7326, 0.01));
    CHECK_THAT(magnetizingCurrentMinimum, Catch::Matchers::WithinRel(0.4347, 0.02));

    settings.reset();
}

// ABT #1002: a moulded body pressed from more than one powder takes each region's own mu. The
// brackets that pin the physics: listing one grade three times must reproduce the single-grade
// model exactly (the regions sum to the piece), a stiffer post must raise the inductance but
// never past the all-stiff body, a plastic-bobbin post must lower it and land on the analytic
// series sum, and a DC bias must pull it down through each region's own knee.
TEST_CASE("Test_Molded_Per_Region_Inductance", "[physical-model][magnetizing-inductance][molded][abt-1002]") {
    settings.reset();
    clear_databases();
    double numberTurns = 12;
    json shapeJson = {
        {"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "molded"},
        {"aliases", json::array()}, {"name", "MAPI-like 4020"},
        {"dimensions", {
            {"A", {{"nominal", 0.0041}}}, {"B", {{"nominal", 0.0021}}}, {"C", {{"nominal", 0.0041}}},
            {"D", {{"nominal", 0.0014}}}, {"E", {{"nominal", 0.0030}}}, {"F", {{"nominal", 0.0012}}}}}
    };
    json windingData = json::parse(
        R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary",
            "numberParallels": 1, "numberTurns": 12, "wire": "Dummy"}]})");
    OpenMagnetics::Coil winding(windingData);
    MagnetizingInductance magnetizingInductanceModel("ZHANG");

    auto inductanceOf = [&](json material, std::string* methodUsed = nullptr) {
        json coreJson;
        coreJson["functionalDescription"] = {
            {"type", "closedShape"}, {"material", material}, {"shape", shapeJson},
            {"gapping", json::array()}, {"numberStacks", 1}};
        Core core(coreJson);
        core.process_data();
        core.process_gap();
        auto output = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(core, winding, nullptr);
        if (methodUsed) {
            *methodUsed = output.get_method_used();
        }
        return output.get_magnetizing_inductance().get_nominal().value();
    };

    std::string bareMethod;
    std::string listMethod;
    double bare26 = inductanceOf("Kool Mµ 26", &bareMethod);
    double listed26 = inductanceOf(json::array({"Kool Mµ 26", "Kool Mµ 26", "Kool Mµ 26"}), &listMethod);
    double paired26 = inductanceOf(json::array({"Kool Mµ 26", "Kool Mµ 26"}));
    UNSCOPED_INFO("bare " << bare26 * 1e6 << " uH (" << bareMethod << "), listed x3 " << listed26 * 1e6
                          << " uH (" << listMethod << "), listed x2 " << paired26 * 1e6 << " uH");
    CHECK(bareMethod != "MoldedPerRegionReluctance");
    CHECK(listMethod == "MoldedPerRegionReluctance");
    CHECK_THAT(listed26, Catch::Matchers::WithinRel(bare26, 1e-6));
    CHECK_THAT(paired26, Catch::Matchers::WithinRel(bare26, 1e-6));

    double bare60 = inductanceOf("Kool Mµ 60");
    double stiffPost = inductanceOf(json::array({"Kool Mµ 60", "Kool Mµ 26"}));
    double stiffBase = inductanceOf(json::array({"Kool Mµ 26", "Kool Mµ 26", "Kool Mµ 60"}));
    CHECK(stiffPost > bare26);
    CHECK(stiffPost < bare60);
    CHECK(stiffBase > bare26);
    CHECK(stiffBase < bare60);

    // A coil on a plastic bobbin: the post is air. Analytic series sum from the piece's own
    // region constants and the grade's zero-bias permeability.
    double bobbinPost = inductanceOf(json::array({"air", "Kool Mµ 26", "Kool Mµ 26"}));
    CHECK(bobbinPost < bare26);
    json coreJson;
    coreJson["functionalDescription"] = {
        {"type", "closedShape"}, {"material", "Kool Mµ 26"}, {"shape", shapeJson},
        {"gapping", json::array()}, {"numberStacks", 1}};
    Core core(coreJson);
    auto regions = CorePiece::factory(core.resolve_shape())->get_region_shape_constants().value();
    // Evaluated at the reference frequency the standard path uses when no operating point is given.
    double referenceFrequency = Defaults().coreAdviserFrequencyReference;
    double permeability26 = InitialPermeability::get_initial_permeability(core.resolve_material(), 25.0, std::nullopt, referenceFrequency);
    double vacuumPermeability = Constants().vacuumPermeability;
    double expectedReluctance = regions[0].c1 / vacuumPermeability +
                                (regions[1].c1 + regions[2].c1) / (vacuumPermeability * permeability26);
    CHECK_THAT(bobbinPost, Catch::Matchers::WithinRel(pow(numberTurns, 2) / expectedReluctance, 1e-6));

    // Four entries do not describe this body.
    CHECK_THROWS(inductanceOf(json::array({"Kool Mµ 26", "Kool Mµ 26", "Kool Mµ 26", "Kool Mµ 26"})));

    // DC bias pulls each region down its own mu(H); the post, with the smallest section, first.
    json biasedCoreJson;
    biasedCoreJson["functionalDescription"] = {
        {"type", "closedShape"}, {"material", json::array({"Kool Mµ 60", "Kool Mµ 26"})}, {"shape", shapeJson},
        {"gapping", json::array()}, {"numberStacks", 1}};
    Core biasedCore(biasedCoreJson);
    biasedCore.process_data();
    biasedCore.process_gap();
    double unbiased = MagnetizingInductance::calculate_molded_magnetizing_inductance(biasedCore, numberTurns, 25, referenceFrequency);
    double biased = MagnetizingInductance::calculate_molded_magnetizing_inductance(biasedCore, numberTurns, 25, referenceFrequency, 1.0);
    // 1 A and 2 A put B = I L / (N A_post) = 0.13 T and 0.22 T in the Kool Mu 60 post (B_sat about 1.0 T),
    // where the flux still rises with the current. It was 20 A and 60 A: from about 3 A on, this per-region
    // iteration's flux FALLS with current (ABT #1765), and from about 10 A its iterate reaches 71 kA/m,
    // where the fit gives mu_r 0.93, below the vacuum's, which the permeability model now refuses as a
    // saturation (ABT #1740).
    double moreBiased = MagnetizingInductance::calculate_molded_magnetizing_inductance(biasedCore, numberTurns, 25, referenceFrequency, 2.0);
    UNSCOPED_INFO("unbiased " << unbiased * 1e6 << " uH, 1 A " << biased * 1e6 << " uH, 2 A " << moreBiased * 1e6 << " uH");
    CHECK(biased < unbiased);
    CHECK(moreBiased < biased);
    CHECK_THAT(unbiased, Catch::Matchers::WithinRel(stiffPost, 1e-6));

    // Through the public entry point WITH an operating point, the way every consumer calls it:
    // a raw current waveform (no processed data, no harmonics) whose DC component is the bias.
    // A symmetric ripple reproduces the unbiased figure; a 1 A offset reproduces the biased one.
    auto operatingPointWithOffset = [&](double offset) {
        json operatingPointJson = json::parse(R"({"conditions": {"ambientTemperature": 25.0},
            "excitationsPerWinding": [{"name": "Primary", "frequency": 100000,
              "current": {"waveform": {"data": [-0.05, 0.05, -0.05], "time": [0, 5e-6, 1e-5]}},
              "voltage": {"waveform": {"data": [-1, 1, 1, -1, -1], "time": [0, 0, 5e-6, 5e-6, 1e-5]}}}]})");
        for (auto& sample : operatingPointJson["excitationsPerWinding"][0]["current"]["waveform"]["data"]) {
            sample = sample.get<double>() + offset;
        }
        return OperatingPoint(operatingPointJson);
    };
    auto rippleOnly = operatingPointWithOffset(0);
    auto withOffset = operatingPointWithOffset(1);
    auto rippleOutput = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(biasedCore, winding, &rippleOnly);
    auto offsetOutput = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(biasedCore, winding, &withOffset);
    CHECK(rippleOutput.get_method_used() == "MoldedPerRegionReluctance");
    CHECK_THAT(rippleOutput.get_magnetizing_inductance().get_nominal().value(), Catch::Matchers::WithinRel(unbiased, 1e-3));
    CHECK_THAT(offsetOutput.get_magnetizing_inductance().get_nominal().value(), Catch::Matchers::WithinRel(biased, 1e-3));
    settings.reset();
}

// ABT #1470: a DC bias that rolls a powder core's permeability off makes L grow far slower than N^2.
// The turns solve assumed N^2, stopped after a fixed Newton budget and returned 118,109 turns for
// 22 uH on a sendust T 2.03/1.27/0.64 at 10 A dc, where the model's inductance was 7.1 uH (below
// the 17.6 uH minimum). A consumer wound that count and took the machine down (#1471). It must throw.
TEST_CASE("Test_Number_Turns_From_Inductance_Unreachable_Under_Dc_Bias_Throws", "[physical-model][magnetizing-inductance][abt-1470]") {
    const double frequency = 100e3;
    const int samples = 256;
    json time = json::array(), data = json::array();
    for (int k = 0; k <= samples; ++k) {   // buck-like triangle: 10 A dc, 4 A peak to peak
        time.push_back(k / (samples * frequency));
        data.push_back(10.0 + (k < samples / 2 ? -2.0 + 4.0 * k / (samples / 2) : 2.0 - 4.0 * (k - samples / 2) / (samples / 2)));
    }
    OpenMagnetics::Inputs inputs(json{
        {"designRequirements", {{"magnetizingInductance", {{"nominal", 22e-6}, {"minimum", 17.6e-6}, {"maximum", 26.4e-6}}},
                                {"turnsRatios", json::array()}, {"isolationSides", {"primary"}}, {"topology", "buckConverter"}}},
        {"operatingPoints", {{{"name", "op 0"}, {"conditions", {{"ambientTemperature", 25}}},
                              {"excitationsPerWinding", {{{"name", "Primary"}, {"frequency", frequency},
                                                          {"current", {{"waveform", {{"time", time}, {"data", data}}}}}}}}}}}});
    for (auto& operatingPoint : inputs.get_mutable_operating_points()) {
        for (auto& excitation : operatingPoint.get_mutable_excitations_per_winding()) {
            auto current = excitation.get_current().value();
            current.set_processed(OpenMagnetics::Inputs::calculate_basic_processed_data(current.get_waveform().value()));
            excitation.set_current(current);
        }
    }
    // ABT #1496: magnetic_autocomplete now throws when the winding does not fit its window, and a
    // 1.25 mm wire does not fit the 0.635 mm ring of T 2.03/1.27/0.64; the tiny toroids are wound
    // with a wire they hold, so the inductance solve is what is under test (Alf, 2026-10-11).
    auto magneticFor = [](const std::string& shape, const std::string& wireName) {
        return magnetic_autocomplete(OpenMagnetics::Magnetic(json{
            {"core", {{"name", shape}, {"functionalDescription", {{"type", "toroidal"}, {"shape", shape}, {"material", "75-Series 26"},
                                                                  {"gapping", json::array()}, {"numberStacks", 1}}}}},
            {"coil", {{"bobbin", "Basic"}, {"functionalDescription", {{{"name", "Primary"}, {"numberTurns", 1}, {"numberParallels", 1},
                                                                        {"isolationSide", "primary"}, {"wire", wireName}}}}}}}));
    };
    MagnetizingInductance magnetizingInductance(ReluctanceModels::ZHANG);

    for (const std::string shape : {"T 2.03/1.27/0.64", "T 3.94/2.24/1.8"}) {
        INFO(shape);
        auto magnetic = magneticFor(shape, "Round 0.2 - Grade 1");
        auto start = std::chrono::steady_clock::now();
        CHECK_THROWS_AS(magnetizingInductance.calculate_number_turns_from_gapping_and_inductance(magnetic.get_core(), magnetic.get_coil(), &inputs),
                        CalculationException);
        CHECK_THROWS_AS(magnetizingInductance.calculate_number_turns_from_gapping_and_inductance(magnetic.get_core(), magnetic.get_coil(), &inputs, DimensionalValues::MINIMUM),
                        CalculationException);
        CHECK(std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() < 30);
    }

    // Control: a core that does reach the target keeps an answer, and that answer meets it.
    auto magnetic = magneticFor("T 21/12/14", "Round 1.25 - Grade 1");
    int numberTurns = magnetizingInductance.calculate_number_turns_from_gapping_and_inductance(magnetic.get_core(), magnetic.get_coil(), &inputs);
    auto coil = magnetic.get_coil();
    coil.get_mutable_functional_description()[0].set_number_turns(numberTurns);
    auto operatingPoint = inputs.get_operating_point(0);
    double inductance = resolve_dimensional_values(
        magnetizingInductance.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), coil, &operatingPoint).get_magnetizing_inductance());
    INFO(numberTurns << " turns, " << inductance * 1e6 << " uH");
    CHECK(inductance >= 17.6e-6);
    CHECK(inductance <= 26.4e-6);
}

TEST_CASE("Test_Gap_From_Inductance_Refuses_A_Bias_That_Collapses_The_Permeability", "[physical-model][magnetizing-inductance][gapping][abt-1855]") {
    // ABT #1855 (Alf, 2026-10-14): a design operation refuses a solution the saturation check flags. 3C95's
    // permeability model falls to the vacuum's at a field strength where its magnetisation curve still
    // carries less than B_sat (criterion 2). Sizing a gap for a DC bias there used to throw from the
    // permeability (main's SaturationException); with mu_r = 1 taken past the knee it would return the
    // gap of an air core instead.
    settings.reset();
    clear_databases();
    const double ambientTemperature = 25;
    const double frequency = 100000;
    const double numberTurns = 20;
    const double desiredMagnetizingInductance = 100e-6;
    const std::string coreMaterial = "3C95";
    Core core;
    OpenMagnetics::Coil coil;
    OpenMagnetics::Inputs inputs;
    prepare_test_parameters(0, ambientTemperature, frequency, numberTurns, desiredMagnetizingInductance, {}, "E 42/21/20",
                            coreMaterial, core, coil, inputs);
    const auto material = core.resolve_material();
    std::optional<double> collapseField;
    for (double field = 10; field < 1e6; field *= 1.05) {
        if (InitialPermeability::get_initial_permeability_model_value(material, ambientTemperature, field, frequency) <= 1) {
            collapseField = field;
            break;
        }
    }
    REQUIRE(collapseField);
    const double collapseFluxDensity = InitialPermeability::get_flux_density_for_magnetic_field_dc_bias(material, collapseField.value(), ambientTemperature, frequency);
    const double saturationFluxDensity = Core::get_magnetic_flux_density_saturation(material, ambientTemperature, false);
    INFO("3C95 permeability model collapses at H = " << collapseField.value() << " A/m, B = " << collapseFluxDensity
         << " T on its curve, B_sat " << saturationFluxDensity << " T");
    REQUIRE(collapseFluxDensity < 0.97 * saturationFluxDensity);
    const double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    // B = L I / (N A_e) at the solution; a 0.1 A ripple keeps the peak within B_sat.
    auto inputsWithBiasFluxDensity = [&](double biasFluxDensity) {
        OpenMagnetics::Coil unused;
        Core unusedCore;
        OpenMagnetics::Inputs biasedInputs;
        prepare_test_parameters(biasFluxDensity * numberTurns * effectiveArea / desiredMagnetizingInductance, ambientTemperature, frequency,
                                numberTurns, desiredMagnetizingInductance, {}, "E 42/21/20", coreMaterial, unusedCore, unused, biasedInputs, 0.1);
        return biasedInputs;
    };
    MagnetizingInductance magnetizingInductance("ZHANG");

    auto collapsing = inputsWithBiasFluxDensity(0.5 * (collapseFluxDensity + 0.97 * saturationFluxDensity));
    REQUIRE_THROWS_AS(magnetizingInductance.calculate_gapping_from_number_turns_and_inductance(core, coil, &collapsing, GappingType::GROUND),
                      SaturationException);
    REQUIRE_THROWS_WITH(magnetizingInductance.calculate_gapping_from_number_turns_and_inductance(core, coil, &collapsing, GappingType::GROUND),
                        Catch::Matchers::ContainsSubstring("permeability model"));

    // Control: half that bias keeps a real gap.
    auto moderate = inputsWithBiasFluxDensity(0.5 * collapseFluxDensity);
    std::vector<CoreGap> gapping;
    REQUIRE_NOTHROW(gapping = magnetizingInductance.calculate_gapping_from_number_turns_and_inductance(core, coil, &moderate, GappingType::GROUND));
    REQUIRE(!gapping.empty());
    CHECK(gapping[0].get_length() > 0);
}

TEST_CASE("Test_Dc_Bias_Is_The_Time_Average_Not_The_Midpoint", "[physical-model][magnetizing-inductance][dc-bias][abt-1223]") {
    // A voltage with zero volt-second average carries no DC, whatever its shape. With a dead time
    // (+V for T/4, -V for T/4, 0 for T/2) the flux is a triangle followed by a flat stretch at its
    // minimum: its time average sits a quarter of the swing above the minimum, its midpoint,
    // (max + min) / 2, half the swing above it. With the DC taken as the average the midpoint is
    // swing / 4 away from zero; read as a bias, it would drive the permeability of this ungapped
    // Nanoperm toroid down the curve. (A plain triangle cannot test this: whatever its duty, its
    // midpoint IS its average, ABT #1467.) The inductance must equal the symmetric square wave's and
    // the unbiased small-signal value.
    auto core = Core(json::parse(R"({"functionalDescription": {"gapping": [], "material": "Nanoperm 80000", "numberStacks": 1, "shape": {"magneticCircuit": "closed", "type": "custom", "family": "t", "aliases": [], "name": "T 63/50/30", "dimensions": {"A": {"nominal": 0.063}, "B": {"nominal": 0.050}, "C": {"nominal": 0.030}}}, "type": "toroidal"}, "manufacturerInfo": null, "name": "My Core"})"));
    auto coil = OpenMagnetics::Coil(json::parse(R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary", "numberParallels": 1, "numberTurns": 1, "wire": "Dummy"}]})"));
    const double period = 1e-5;
    auto operatingPointWith = [period](std::vector<double> data, std::vector<double> time) {
        json excitationJson;
        excitationJson["frequency"] = 1 / period;
        excitationJson["name"] = "Primary";
        excitationJson["voltage"]["waveform"]["data"] = data;
        excitationJson["voltage"]["waveform"]["time"] = time;
        json operatingPointJson;
        operatingPointJson["conditions"]["ambientTemperature"] = 25.0;
        operatingPointJson["excitationsPerWinding"] = json::array({excitationJson});
        return OperatingPoint(operatingPointJson);
    };
    MagnetizingInductance magnetizingInductance(ReluctanceModels::ZHANG);

    // Same volt-seconds per excursion in both (68.85 V for T/4 = 34.425 V for T/2), so the same flux
    // swing, about 1.1 T peak to peak: below Nanoperm's saturation.
    const double voltage = 68.85;
    auto deadTime = operatingPointWith({voltage, voltage, -voltage, -voltage, 0, 0},
                                       {0, period / 4, period / 4, period / 2, period / 2, period});
    auto symmetric = operatingPointWith({voltage / 2, voltage / 2, -voltage / 2, -voltage / 2, voltage / 2},
                                        {0, period / 2, period / 2, period, period});
    auto [deadTimeOutput, deadTimeFluxDensity] = magnetizingInductance.calculate_inductance_and_magnetic_flux_density(core, coil, &deadTime);
    auto [symmetricOutput, symmetricFluxDensity] = magnetizingInductance.calculate_inductance_and_magnetic_flux_density(core, coil, &symmetric);
    double deadTimeInductance = resolve_dimensional_values(deadTimeOutput.get_magnetizing_inductance());
    double symmetricInductance = resolve_dimensional_values(symmetricOutput.get_magnetizing_inductance());

    auto fluxData = deadTimeFluxDensity.get_waveform().value().get_data();
    double swing = *std::max_element(fluxData.begin(), fluxData.end()) - *std::min_element(fluxData.begin(), fluxData.end());
    double midpoint = deadTimeFluxDensity.get_processed().value().get_offset();
    double average = OpenMagnetics::Inputs::calculate_waveform_average(deadTimeFluxDensity.get_waveform().value());
    INFO("dead-time flux: swing " << swing << " T, midpoint " << midpoint << " T, average " << average << " T; L " << deadTimeInductance * 1e6 << " uH vs " << symmetricInductance * 1e6 << " uH");
    // The case is meaningful only if the midpoint is far from the average: a quarter of the swing.
    REQUIRE(swing > 0.5);
    // (The shape alone: midpoint - average does not depend on the DC the model assigns.)
    REQUIRE_THAT(midpoint - average, Catch::Matchers::WithinAbs(swing / 4, 0.02 * swing));

    double unbiasedPermeability = InitialPermeability().get_initial_permeability(core.resolve_material(), 25.0, std::nullopt, 1e5);
    double unbiasedInductance = 1.0 / ReluctanceModel::factory(ReluctanceModels::ZHANG)->get_core_reluctance(core, unbiasedPermeability).get_core_reluctance();
    INFO("unbiased L " << unbiasedInductance * 1e6 << " uH");
    CHECK_THAT(deadTimeInductance, Catch::Matchers::WithinRel(symmetricInductance, 1e-3));
    CHECK_THAT(deadTimeInductance, Catch::Matchers::WithinRel(unbiasedInductance, 1e-2));
    // Zero net volt-seconds: no DC flux.
    CHECK_THAT(average, Catch::Matchers::WithinAbs(0, 0.01 * swing));
}

TEST_CASE("Test_Volt_Seconds_That_Saturate_The_Core_Are_Refused_Before_The_Dc_Bias_Runs_Away", "[physical-model][magnetizing-inductance][dc-bias][abt-1740]") {
    // The 3C96 Icepak flyback fixture (E 25/13/7, A_e 51.8 mm2, 24:13 on a 0.8 mm gap) used to carry a
    // second operating point of 911 uVs, a 0.73 T swing on a material that saturates at about 0.5 T. The
    // DC-bias fixed point ran away on it (H 224 -> 3.6e17 A/m, mu 309 -> 8e-13) and reported "biased
    // permeability at mu = 0.000000 is 0.000000". That point is gone from the fixture; the same swing is
    // reproduced here from its first operating point (595 uVs, 0.48 T on 24 turns) with 16 turns:
    // 0.48 * 24 / 16 = 0.72 T. It must say the core saturates, naming the operating point and winding.
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "concentric_flyback_rectangular_column.json");
    auto mas = OpenMagneticsTesting::mas_loader(path.string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = mas.get_inputs();
    MagnetizingInductance magnetizingInductance(ReluctanceModels::ZHANG);
    REQUIRE(inputs.get_operating_point(0).get_name().value() == "Min. input volt. (simulated)");
    auto coilWithTurns = [&](int numberTurns) {
        auto coil = magnetic.get_coil();
        coil.get_mutable_functional_description()[0].set_number_turns(numberTurns);
        return coil;
    };

    // Until ABT #1855 both were refused (SaturationException naming the operating point and winding).
    // Since then (Alf approved this rewrite 2026-10-14, ABT #1855) the model gives its numbers for a
    // saturated core, the DC-bias fixed point still settles (it no longer runs away), and the saturation
    // is the verdict of the shared check analyze_results makes, ResultsAnalyzer::check_saturation.
    auto material = magnetic.get_core().resolve_material();
    double temperature = inputs.get_operating_point(0).get_conditions().get_ambient_temperature();
    auto swingOf = [](const SignalDescriptor& fluxDensity) {
        auto data = fluxDensity.get_waveform().value().get_data();
        return *std::max_element(data.begin(), data.end()) - *std::min_element(data.begin(), data.end());
    };

    auto real = inputs.get_operating_point(0);
    auto realFluxDensity = magnetizingInductance.calculate_inductance_and_magnetic_flux_density(magnetic.get_core(), magnetic.get_coil(), &real).second;
    double realSwing = swingOf(realFluxDensity);

    auto saturating = inputs.get_operating_point(0);
    auto start = std::chrono::steady_clock::now();
    std::pair<MagnetizingInductanceOutput, SignalDescriptor> saturatingResult;
    REQUIRE_NOTHROW(saturatingResult = magnetizingInductance.calculate_inductance_and_magnetic_flux_density(magnetic.get_core(), coilWithTurns(16), &saturating));
    CHECK(std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() < 30);
    double saturatingInductance = resolve_dimensional_values(saturatingResult.first.get_magnetizing_inductance());
    double saturatingSwing = swingOf(saturatingResult.second);
    double saturatingPeak = ResultsAnalyzer::get_peak_flux_density(saturatingResult.second);
    auto saturatingVerdict = ResultsAnalyzer::check_saturation(material, saturatingPeak, temperature);
    INFO("16 turns: L " << saturatingInductance * 1e6 << " uH, swing " << saturatingSwing << " T, B peak " << saturatingPeak
         << " T, B_sat " << saturatingVerdict.saturationFluxDensity << " T; 24 turns: swing " << realSwing << " T");
    CHECK(std::isfinite(saturatingInductance));
    CHECK(saturatingInductance > 0);
    // Faraday: the same volt-seconds on 16 turns instead of 24.
    CHECK_THAT(saturatingSwing, Catch::Matchers::WithinRel(realSwing * 24 / 16, 0.01));
    CHECK(saturatingVerdict.saturates);
    CHECK_THAT(saturatingVerdict.message, Catch::Matchers::ContainsSubstring("The core saturates"));

    // 8 turns triple the swing to about 1.4 T peak-to-peak: more than twice B_sat, which no bias can hold.
    auto tooFewTurns = inputs.get_operating_point(0);
    std::pair<MagnetizingInductanceOutput, SignalDescriptor> tooFewTurnsResult;
    REQUIRE_NOTHROW(tooFewTurnsResult = magnetizingInductance.calculate_inductance_and_magnetic_flux_density(magnetic.get_core(), coilWithTurns(8), &tooFewTurns));
    double tooFewTurnsSwing = swingOf(tooFewTurnsResult.second);
    INFO("8 turns: swing " << tooFewTurnsSwing << " T");
    CHECK_THAT(tooFewTurnsSwing, Catch::Matchers::WithinRel(realSwing * 3, 0.01));
    CHECK(tooFewTurnsSwing / 2 > saturatingVerdict.saturationFluxDensity);
    CHECK(ResultsAnalyzer::check_saturation(material, ResultsAnalyzer::get_peak_flux_density(tooFewTurnsResult.second), temperature).saturates);

    // Control: the real coil (24 turns, 0.48 T swing) has a DC point and keeps its answer, the ~64 uH
    // the fixture's design requirements now state.
    auto converging = inputs.get_operating_point(0);
    auto [output, fluxDensity] = magnetizingInductance.calculate_inductance_and_magnetic_flux_density(magnetic.get_core(), magnetic.get_coil(), &converging);
    double inductance = resolve_dimensional_values(output.get_magnetizing_inductance());
    INFO("L " << inductance * 1e6 << " uH, B peak " << fluxDensity.get_processed().value().get_peak().value() << " T");
    CHECK(inductance >= 60e-6);
    CHECK(inductance <= 67e-6);
}

TEST_CASE("Test_Magnetic_Flux_Density_Dcm_Flyback_Keeps_The_Dead_Time", "[physical-model][magnetizing-inductance][abt-1467]") {
    // ABT #1467: the core-loss path routes a flyback-labelled multi-winding operating point to the
    // multiport-inductor branch, which rebuilt the magnetizing current as a triangle rising over D
    // and falling over all of 1 - D. A DCM flyback's flux is flat through the dead time (0 V on
    // every winding for the last 42 of 128 samples); Faraday on the primary voltage keeps it. The
    // fixture is the #1463 flyback with its reset balanced (-235 * 32 / 54 V instead of -137.2 V):
    // the mean removed from an unbalanced voltage would tilt the dead time.
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "Test_Magnetizing_Current_Dcm_Flyback_Balanced.json");
    std::ifstream in(path);
    json masJson = json::parse(in);
    OpenMagnetics::Magnetic magnetic(masJson["magnetic"]);
    for (auto& excitation : masJson["inputs"]["operatingPoints"][0]["excitationsPerWinding"]) {
        excitation.erase("magnetizingCurrent");
    }
    OperatingPoint operatingPoint(masJson["inputs"]["operatingPoints"][0]);
    MagnetizingInductance magnetizingInductance;
    auto [output, magneticFluxDensity] = magnetizingInductance.calculate_inductance_and_magnetic_flux_density(magnetic.get_core(), magnetic.get_coil(), &operatingPoint);
    // Stored compressed (corner points only): sample it over the period.
    auto sampledFluxDensity = OpenMagnetics::Inputs::calculate_sampled_waveform(magneticFluxDensity.get_waveform().value(), 100000);
    const auto data = sampledFluxDensity.get_data();
    REQUIRE(data.size() >= 128);
    double maximum = *std::max_element(data.begin(), data.end());
    double minimum = *std::min_element(data.begin(), data.end());
    double deadTimeLevel = data.back();
    size_t flatSamples = 0;
    for (auto value : data) {
        if (std::fabs(value - deadTimeLevel) <= 1e-6 * (maximum - minimum)) {
            flatSamples++;
        }
    }
    double flatFraction = static_cast<double>(flatSamples) / data.size();
    // Faraday: B swing = V t_on / (N_p A_e).
    double effectiveArea = magnetic.get_core().get_processed_description()->get_effective_parameters().get_effective_area();
    double faradaySwing = 235.0 * 2.5e-6 / (42 * effectiveArea);
    INFO("B max " << maximum << " T, min " << minimum << " T, flat fraction " << flatFraction << ", Faraday swing " << faradaySwing << " T");
    CHECK(flatFraction >= 0.3);
    CHECK_THAT(maximum - minimum, Catch::Matchers::WithinRel(faradaySwing, 0.03));
}

TEST_CASE("Test_Multiport_Flyback_Without_Voltage_Is_The_Ampere_Turn_Sum", "[physical-model][magnetizing-inductance][abt-1467]") {
    // ABT #1467: a flyback-labelled operating point with no primary voltage. Its magnetizing current
    // is the windings' ampere-turn sum referred to the primary (secondary signed by its isolation
    // side), exact while either winding conducts: it keeps the DCM dead time (both currents zero)
    // where the FLYBACK label's triangle fell over all of 1 - D.
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "Test_Winding_Losses_Abt1463_Flyback.json");
    std::ifstream in(path);
    json masJson = json::parse(in);
    for (auto& excitation : masJson["inputs"]["operatingPoints"][0]["excitationsPerWinding"]) {
        excitation.erase("magnetizingCurrent");
        excitation.erase("voltage");
    }
    OperatingPoint operatingPoint(masJson["inputs"]["operatingPoints"][0]);
    auto magnetizingCurrent = OpenMagnetics::Inputs::get_multiport_inductor_magnetizing_current(
        operatingPoint, 250e-6, {7.0}, std::vector<IsolationSide>{IsolationSide::PRIMARY, IsolationSide::SECONDARY});
    auto waveform = magnetizingCurrent.get_waveform().value();
    const auto& data = waveform.get_data();
    double maximum = *std::max_element(data.begin(), data.end());
    double minimum = *std::min_element(data.begin(), data.end());
    size_t zeroSamples = 0;
    for (auto value : data) {
        if (std::fabs(value) <= 1e-9 * maximum) {
            zeroSamples++;
        }
    }
    INFO("i_m max " << maximum << " A, min " << minimum << " A, " << zeroSamples << " of " << data.size() << " samples at zero");
    // Peak: the primary's 128-sample peak, 2.277 A; the secondary's 16.16 A / 7 = 2.309 A at commutation.
    CHECK_THAT(maximum, Catch::Matchers::WithinRel(2.31, 0.02));
    CHECK(minimum > -1e-9 * maximum);
    // Dead time: the last 42 of 128 samples carry no current on either winding.
    CHECK(static_cast<double>(zeroSamples) / data.size() >= 40.0 / 128);
}

TEST_CASE("Test_Forward_Transformer_Load_Current_Is_Not_Core_Bias", "[physical-model][magnetizing-inductance][abt-1467][abt-1814]") {
    // ABT #1814, fixed by ABT #1467: MAS example 04 (forward transformer, E 32/16 N87, 3 windings).
    // The primary's load current used to be taken as the core's DC bias: B offset 0.4055 T and an
    // inductance of 34.304 uH at the operating point, 5 % of the unbiased 673.212 uH. The magnetizing
    // current is Faraday on the winding voltage, so the load ampere-turns, which the secondary
    // cancels, do not bias the core: |B offset| is about 3.5 mT and the inductance is the unbiased one.
    settings.reset();
    clear_databases();
    std::ifstream file(std::string(__FILE__).substr(0, std::string(__FILE__).rfind('/'))
                       + "/../MAS/examples/04_forward_xfmr_e3216_n87.json");
    REQUIRE(file.is_open());
    json masJson = json::parse(file);
    OpenMagnetics::Magnetic magnetic(masJson["magnetic"]);
    magnetic = magnetic_autocomplete(magnetic);
    for (auto& operatingPointJson : masJson["inputs"]["operatingPoints"]) {
        for (auto& excitationJson : operatingPointJson["excitationsPerWinding"]) {
            excitationJson.erase("magnetizingCurrent");
        }
    }
    OpenMagnetics::Inputs inputs(masJson["inputs"], false);
    REQUIRE(inputs.get_operating_points().size() >= 1);
    OperatingPoint operatingPoint = inputs.get_operating_points()[0];
    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 3);

    MagnetizingInductance model;
    double unbiasedInductance = resolve_dimensional_values(
        model.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil()).get_magnetizing_inductance());
    auto [output, magneticFluxDensity] = model.calculate_inductance_and_magnetic_flux_density(magnetic.get_core(), magnetic.get_coil(), &operatingPoint);
    double inductanceAtOperatingPoint = resolve_dimensional_values(output.get_magnetizing_inductance());

    CHECK_THAT(unbiasedInductance, Catch::Matchers::WithinRel(673.212e-6, 0.01));
    CHECK_THAT(inductanceAtOperatingPoint, Catch::Matchers::WithinRel(unbiasedInductance, 0.01));
    CHECK(std::abs(magneticFluxDensity.get_processed()->get_offset()) < 0.01);
}

TEST_CASE("Test_Forward_And_Multi_Operating_Point_Transformers_Load_Current_Is_Not_Core_Bias", "[physical-model][magnetizing-inductance][abt-1814]") {
    // ABT #1814, the two other MAS examples the collapse hit (only 04 was guarded): the primary's
    // reflected load current was taken as the core's DC bias. 22 (forward, E 42/21/20 3C95, 2
    // windings) fell from 3438 uH to 359.7 uH with a B offset of 0.377 T; 19 (2 operating points,
    // ETD 39/20 PC95) fell at its first operating point from 1583 uH to 954 uH with 0.302 T. Their
    // windings' net ampere-turns average zero, so the B offset is a few mT at every operating point
    // and the inductance stays near the unbiased one (19's sits 1.7 % above it: the permeability at
    // its small flux swing and temperature).
    struct Example { std::string file; double unbiasedInductance; double tolerance; };
    for (const auto& example : {Example{"22_margin_tape_forward_e4218_3c95", 3438.18e-6, 0.01},
                                Example{"19_multi_op_xfmr_etd3920_pc95", 1583.49e-6, 0.03}}) {
        settings.reset();
        clear_databases();
        std::ifstream file(std::string(__FILE__).substr(0, std::string(__FILE__).rfind('/'))
                           + "/../MAS/examples/" + example.file + ".json");
        REQUIRE(file.is_open());
        json masJson = json::parse(file);
        OpenMagnetics::Magnetic magnetic(masJson["magnetic"]);
        magnetic = magnetic_autocomplete(magnetic);
        for (auto& operatingPointJson : masJson["inputs"]["operatingPoints"]) {
            for (auto& excitationJson : operatingPointJson["excitationsPerWinding"]) {
                excitationJson.erase("magnetizingCurrent");
            }
        }
        OpenMagnetics::Inputs inputs(masJson["inputs"], false);
        REQUIRE(inputs.get_operating_points().size() >= 1);

        MagnetizingInductance model;
        double unbiasedInductance = resolve_dimensional_values(
            model.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil()).get_magnetizing_inductance());
        CHECK_THAT(unbiasedInductance, Catch::Matchers::WithinRel(example.unbiasedInductance, 0.01));
        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
            OperatingPoint operatingPoint = inputs.get_operating_points()[operatingPointIndex];
            REQUIRE(operatingPoint.get_excitations_per_winding().size() >= 2);
            auto [output, magneticFluxDensity] = model.calculate_inductance_and_magnetic_flux_density(magnetic.get_core(), magnetic.get_coil(), &operatingPoint);
            double inductanceAtOperatingPoint = resolve_dimensional_values(output.get_magnetizing_inductance());
            INFO(example.file << " operating point " << operatingPointIndex << ": L " << inductanceAtOperatingPoint * 1e6
                 << " uH (unbiased " << unbiasedInductance * 1e6 << " uH), B offset " << magneticFluxDensity.get_processed()->get_offset() << " T");
            CHECK_THAT(inductanceAtOperatingPoint, Catch::Matchers::WithinRel(unbiasedInductance, example.tolerance));
            CHECK(std::abs(magneticFluxDensity.get_processed()->get_offset()) < 0.01);
        }
    }
}

TEST_CASE("Test_Magnetizing_Current_Throws_When_A_Winding_Has_No_Current", "[physical-model][magnetizing-inductance][abt-1814]") {
    // ABT #1814 (3): the DC of a multi-winding magnetizing current comes from the sum of every
    // winding's ampere-turns. With a current on the primary and none on a secondary that sum
    // cannot be formed; the conduction map used to come back empty and the magnetizing current
    // silently lost its DC. It throws now, naming the winding. An operating point given by its
    // voltages alone (no winding with a current) stays a voltage-only one: no DC from currents.
    settings.reset();
    clear_databases();
    std::ifstream file(std::string(__FILE__).substr(0, std::string(__FILE__).rfind('/'))
                       + "/../MAS/examples/04_forward_xfmr_e3216_n87.json");
    REQUIRE(file.is_open());
    json masJson = json::parse(file);
    OpenMagnetics::Magnetic magnetic(masJson["magnetic"]);
    magnetic = magnetic_autocomplete(magnetic);
    for (auto& excitationJson : masJson["inputs"]["operatingPoints"][0]["excitationsPerWinding"]) {
        excitationJson.erase("magnetizingCurrent");
    }
    json partialJson = masJson["inputs"]["operatingPoints"][0];
    partialJson["excitationsPerWinding"][2].erase("current");
    OperatingPoint partial(partialJson);
    REQUIRE(partial.get_excitations_per_winding()[0].get_current());
    REQUIRE_FALSE(partial.get_excitations_per_winding()[2].get_current());

    CHECK_THROWS_WITH(OpenMagnetics::Inputs::calculate_winding_conduction(partial),
                      Catch::Matchers::ContainsSubstring("winding 2 has no current"));
    CHECK_THROWS_WITH(OpenMagnetics::Inputs::calculate_magnetizing_current(partial, 0, 673e-6, magnetic.get_coil().get_turns_ratios(), false),
                      Catch::Matchers::ContainsSubstring("winding 2 has no current"));
    CHECK_THROWS_WITH(OpenMagnetics::Inputs::include_dc_offset_into_magnetizing_current(partial, magnetic.get_coil().get_turns_ratios()),
                      Catch::Matchers::ContainsSubstring("winding 2 has no current"));

    json voltageOnlyJson = masJson["inputs"]["operatingPoints"][0];
    for (auto& excitationJson : voltageOnlyJson["excitationsPerWinding"]) {
        excitationJson.erase("current");
    }
    OperatingPoint voltageOnly(voltageOnlyJson);
    CHECK_FALSE(OpenMagnetics::Inputs::calculate_winding_conduction(voltageOnly).has_value());
    CHECK_NOTHROW(OpenMagnetics::Inputs::calculate_magnetizing_current(voltageOnly, 0, 673e-6, magnetic.get_coil().get_turns_ratios(), false));
}

TEST_CASE("Test_Gapping_Reads_The_Magnetizing_Current_Of_The_Target_Inductance", "[physical-model][magnetizing-inductance][gapping][abt-1902]") {
    // ABT #1902: a forward transformer carried in a magnetizing current a converter model had derived from the
    // winding voltage for its minimum inductance; the gap solve read that current against the reluctance of the
    // (eight times larger) target and found a B peak eight times the volt-seconds' one: "no gap holds it, the
    // core saturates". With a voltage the flux is Faraday's, and the gap solve derives the magnetizing current
    // at the target inductance, as the inductance model does: the carried current does not change the gaps.
    settings.reset();
    clear_databases();
    std::ifstream file(std::string(__FILE__).substr(0, std::string(__FILE__).rfind('/'))
                       + "/../MAS/examples/04_forward_xfmr_e3216_n87.json");
    REQUIRE(file.is_open());
    json masJson = json::parse(file);
    OpenMagnetics::Magnetic magnetic(masJson["magnetic"]);
    magnetic = magnetic_autocomplete(magnetic);
    for (auto& operatingPointJson : masJson["inputs"]["operatingPoints"]) {
        for (auto& excitationJson : operatingPointJson["excitationsPerWinding"]) {
            excitationJson.erase("magnetizingCurrent");
        }
    }
    OpenMagnetics::Inputs inputs(masJson["inputs"], false);
    REQUIRE(inputs.get_operating_points().size() >= 1);
    REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage());

    MagnetizingInductance model;
    double unbiasedInductance = resolve_dimensional_values(
        model.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil()).get_magnetizing_inductance());
    double targetInductance = unbiasedInductance / 2;
    DimensionWithTolerance target;
    target.set_nominal(targetInductance);
    inputs.get_mutable_design_requirements().set_magnetizing_inductance(target);

    auto inputsWithCurrentAt = [&](double inductance) {
        OpenMagnetics::Inputs stamped = inputs;
        auto operatingPoint = stamped.get_operating_point(0);
        auto magnetizingCurrent = OpenMagnetics::Inputs::calculate_magnetizing_current(operatingPoint, 0, inductance,
                                                                                         magnetic.get_coil().get_turns_ratios(), false);
        operatingPoint.get_mutable_excitations_per_winding()[0].set_magnetizing_current(magnetizingCurrent);
        stamped.set_operating_point_by_index(operatingPoint, 0);
        return stamped;
    };
    auto staleInputs = inputsWithCurrentAt(targetInductance / 8);
    auto targetInputs = inputsWithCurrentAt(targetInductance);

    auto gapsFromTarget = model.calculate_gapping_from_number_turns_and_inductance(magnetic.get_core(), magnetic.get_coil(), &targetInputs, GappingType::GROUND, 6);
    std::vector<CoreGap> gapsFromStale;
    REQUIRE_NOTHROW(gapsFromStale = model.calculate_gapping_from_number_turns_and_inductance(magnetic.get_core(), magnetic.get_coil(), &staleInputs, GappingType::GROUND, 6));
    REQUIRE(gapsFromStale.size() == gapsFromTarget.size());
    for (size_t i = 0; i < gapsFromStale.size(); ++i) {
        CHECK(gapsFromStale[i].get_length() == gapsFromTarget[i].get_length());
    }

    Core gappedCore = magnetic.get_core();
    gappedCore.set_processed_description(std::nullopt);
    gappedCore.set_geometrical_description(std::nullopt);
    gappedCore.get_mutable_functional_description().set_gapping(gapsFromStale);
    gappedCore.process_data();
    gappedCore.process_gap();
    double gappedInductance = resolve_dimensional_values(
        model.calculate_inductance_from_number_turns_and_gapping(gappedCore, magnetic.get_coil()).get_magnetizing_inductance());
    CHECK_THAT(gappedInductance, Catch::Matchers::WithinRel(targetInductance, 0.05));
}

// ABT #1872: calculate_gap_from_saturation_constraint grows the gap 15 times by 1.5x; when B
// is still above the target after the last step it used to return the 16th gap, one it had
// never evaluated, as if it met the target. With 0.1 nH at 1 A on an ETD 39 the turns solve
// stays at one turn, so B = I / (R A) is bounded by the ungapped core reluctance plus a
// last tried 1 um x 1.5^14 gap (~4 mT) and never reaches 1 mT: the search must say so.
TEST_CASE("Test_Gap_From_Saturation_Constraint_Throws_When_Not_Converged", "[physical-model][magnetizing-inductance][gapping][abt-1872]") {
    settings.reset();
    clear_databases();
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(100000, 1e-10, 25, WaveformLabel::SINUSOIDAL, 2, 0.5, 0);
    auto core = OpenMagneticsTesting::get_quick_core("ETD 39", OpenMagneticsTesting::get_residual_gap(), 1, "N87");
    MagnetizingInductance model;

    REQUIRE_THROWS_AS(model.calculate_gap_from_saturation_constraint(core, &inputs, 0.001, 1.0), CalculationException);

    // A target the search can reach still returns a gap.
    double gap = 0;
    REQUIRE_NOTHROW(gap = model.calculate_gap_from_saturation_constraint(core, &inputs, 0.3, 1.0));
    CHECK(gap > 0);
}

// ABT #1814 (4), Alf 2026-10-20: "use voltage for anything with more than one winding, and current for inductors".
// Which excitation drives the core is decided by the circuit's winding count, never by who stamped a carried
// magnetizing current (a stamp round-tripped through JSON cannot be told from the caller's own).
namespace {
    double magnetizing_current_peak(const OperatingPoint& operatingPoint) {
        auto magnetizingCurrent = operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current();
        REQUIRE(magnetizingCurrent);
        double peak = 0;
        const auto data = magnetizingCurrent->get_waveform().value().get_data();
        for (double value : data) {
            peak = std::max(peak, std::abs(value));
        }
        return peak;
    }
}

TEST_CASE("Test_Transformer_Derives_Its_Magnetizing_Current_From_The_Voltage", "[physical-model][magnetizing-inductance][abt-1814]") {
    // The forward transformer of MAS example 04 carries, with its voltages, a magnetizing current derived for an
    // eighth of its inductance (as a converter model's minimum-Lm current would be), set in memory or read back
    // from JSON: either way the core is driven by the current the voltage gives at the inductance evaluated.
    settings.reset();
    clear_databases();
    std::ifstream file(std::string(__FILE__).substr(0, std::string(__FILE__).rfind('/'))
                       + "/../MAS/examples/04_forward_xfmr_e3216_n87.json");
    REQUIRE(file.is_open());
    json masJson = json::parse(file);
    OpenMagnetics::Magnetic magnetic(masJson["magnetic"]);
    magnetic = magnetic_autocomplete(magnetic);
    json inputsJson = masJson["inputs"];
    for (auto& operatingPointJson : inputsJson["operatingPoints"]) {
        for (auto& excitationJson : operatingPointJson["excitationsPerWinding"]) {
            excitationJson.erase("magnetizingCurrent");
        }
    }

    MagnetizingInductance model;
    double inductance = resolve_dimensional_values(
        model.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil()).get_magnetizing_inductance());

    OpenMagnetics::Inputs inputs(inputsJson, false);
    auto withoutMagnetizingCurrent = inputs.get_operating_point(0);
    REQUIRE(withoutMagnetizingCurrent.get_excitations_per_winding().size() > 1);
    REQUIRE(withoutMagnetizingCurrent.get_excitations_per_winding()[0].get_voltage());
    auto eightfoldCurrent = OpenMagnetics::Inputs::calculate_magnetizing_current(withoutMagnetizingCurrent, 0, inductance / 8,
                                                                                 magnetic.get_coil().get_turns_ratios(), false);

    auto fromVoltage = withoutMagnetizingCurrent;
    model.calculate_inductance_and_magnetic_flux_density(magnetic.get_core(), magnetic.get_coil(), &fromVoltage);
    double voltagePeak = magnetizing_current_peak(fromVoltage);

    auto carried = withoutMagnetizingCurrent;
    carried.get_mutable_excitations_per_winding()[0].set_magnetizing_current(eightfoldCurrent);
    model.calculate_inductance_and_magnetic_flux_density(magnetic.get_core(), magnetic.get_coil(), &carried);

    json carriedJson = inputsJson;
    carriedJson["operatingPoints"][0]["excitationsPerWinding"][0]["magnetizingCurrent"] = eightfoldCurrent;
    OpenMagnetics::Inputs fromJson(carriedJson, false);
    auto readBack = fromJson.get_operating_point(0);
    REQUIRE(readBack.get_excitations_per_winding()[0].get_magnetizing_current());
    model.calculate_inductance_and_magnetic_flux_density(magnetic.get_core(), magnetic.get_coil(), &readBack);

    INFO("magnetizing current peak: from the voltage " << voltagePeak << " A, carried eightfold " << magnetizing_current_peak(carried)
         << " A, read back from JSON " << magnetizing_current_peak(readBack) << " A");
    CHECK_THAT(magnetizing_current_peak(carried), Catch::Matchers::WithinRel(voltagePeak, 1e-6));
    CHECK_THAT(magnetizing_current_peak(readBack), Catch::Matchers::WithinRel(voltagePeak, 1e-6));
}

TEST_CASE("Test_Inductor_Uses_The_Current_It_Is_Given", "[physical-model][magnetizing-inductance][abt-1814]") {
    // A single-winding inductor with a voltage: its given current drives the core, never the volt-seconds. With a
    // winding current that current is the magnetizing current, whatever magnetizing current is carried; with only
    // a magnetizing current, that one is used as it comes.
    settings.reset();
    clear_databases();
    Core core;
    OpenMagnetics::Coil coil;
    OpenMagnetics::Inputs inputs;
    prepare_test_parameters(0, 25, 100000, 40, 200e-6, OpenMagneticsTesting::get_ground_gap(0.0005), "ETD 29", "3C97",
                            core, coil, inputs, 2);
    auto operatingPoint = inputs.get_operating_point(0);
    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 1);
    auto excitation = operatingPoint.get_excitations_per_winding()[0];
    REQUIRE(excitation.get_voltage());
    REQUIRE(excitation.get_current());
    auto tripledCurrent = OpenMagnetics::Inputs::get_excitation_with_proportional_current(excitation, 3).get_current().value();
    double currentPeak = 0;
    const auto currentData = excitation.get_current()->get_waveform().value().get_data();
    for (double value : currentData) {
        currentPeak = std::max(currentPeak, std::abs(value));
    }
    REQUIRE(currentPeak > 0);

    MagnetizingInductance model;

    auto withCurrent = operatingPoint;
    withCurrent.get_mutable_excitations_per_winding()[0].set_magnetizing_current(tripledCurrent);
    model.calculate_inductance_and_magnetic_flux_density(core, coil, &withCurrent);

    auto onlyMagnetizingCurrent = operatingPoint;
    onlyMagnetizingCurrent.get_mutable_excitations_per_winding()[0].set_current(std::nullopt);
    onlyMagnetizingCurrent.get_mutable_excitations_per_winding()[0].set_magnetizing_current(tripledCurrent);
    model.calculate_inductance_and_magnetic_flux_density(core, coil, &onlyMagnetizingCurrent);

    INFO("current peak " << currentPeak << " A; magnetizing current peak with the winding current "
         << magnetizing_current_peak(withCurrent) << " A, with only a tripled magnetizing current "
         << magnetizing_current_peak(onlyMagnetizingCurrent) << " A");
    CHECK_THAT(magnetizing_current_peak(withCurrent), Catch::Matchers::WithinRel(currentPeak, 0.02));
    CHECK_THAT(magnetizing_current_peak(onlyMagnetizingCurrent), Catch::Matchers::WithinRel(3 * currentPeak, 0.02));
}

namespace {
// ABT #1942 fixtures: the MAS file's magnetic and one operating point, nulls dropped as the app's own loader does.
void strip_nulls_abt1942(json& value) {
    if (value.is_array()) {
        for (auto& element : value) {
            strip_nulls_abt1942(element);
        }
        return;
    }
    if (!value.is_object()) {
        return;
    }
    for (auto it = value.begin(); it != value.end();) {
        if (it->is_null()) {
            it = value.erase(it);
        }
        else {
            strip_nulls_abt1942(*it);
            ++it;
        }
    }
}
json load_abt1942(const std::string& fileName) {
    std::ifstream file(OpenMagneticsTesting::get_test_data_path(std::source_location::current(), fileName));
    REQUIRE(file.good());
    json masJson = json::parse(file);
    strip_nulls_abt1942(masJson);
    return masJson;
}
// The gaps' reluctance (with their fringing) of the default reluctance model: the core's total minus the ungapped core's.
double gap_reluctance_abt1942(const Core& core, double permeability) {
    auto reluctance = ReluctanceModel::factory()->get_core_reluctance(core, permeability);
    return reluctance.get_core_reluctance() - reluctance.get_ungapped_core_reluctance().value();
}
std::pair<double, double> extremes_abt1942(const std::vector<double>& data) {
    return {*std::min_element(data.begin(), data.end()), *std::max_element(data.begin(), data.end())};
}
}

TEST_CASE("Test_Dcm_Imposed_Current_Flux_Density_Rises_From_Zero_On_The_Magnetisation_Curve", "[physical-model][magnetizing-inductance][dc-bias][abt-1942]") {
    // ABT #1942: low_filling_inductor (EL 25/6.3, Fair-Rite 78, 6 turns, residual gaps) carries 0 -> 5.73 A -> 0 every
    // 200 kHz cycle, average 2.25 A. One winding, so its current is the magnetizing current: the MMF is imposed. The load
    // line of the average MMF (H 275 A/m) and the ripple linear about it put B at 0.273 -> 0.537 T, a DC floor and a
    // peak above the magnetisation curve. The core comes back to zero MMF every cycle, so B is the curve's from zero:
    // 0 when the current is, B_LL(N * 5.73 A) = 0.439 T at the peak, on the load line N i = H l_e + B A_e R_gap.
    // (The fixture's voltage, 7.29 V for 2.25 us, implies 2.9 uH and contradicts its current; the current decides.)
    settings.reset();
    auto masJson = load_abt1942("low_filling_inductor.json");
    OpenMagnetics::Magnetic magnetic(masJson["magnetic"]);
    OpenMagnetics::Inputs inputs(masJson["inputs"]);
    auto operatingPoint = inputs.get_operating_points()[0];
    auto core = magnetic.get_core();
    auto material = core.resolve_material();
    double numberTurns = magnetic.get_coil().get_functional_description()[0].get_number_turns();
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    double effectiveLength = core.get_processed_description()->get_effective_parameters().get_effective_length();
    double temperature = operatingPoint.get_conditions().get_ambient_temperature();
    double frequency = operatingPoint.get_excitations_per_winding()[0].get_frequency();

    MagnetizingInductance magnetizingInductance;
    auto [output, fluxDensity] = magnetizingInductance.calculate_inductance_and_magnetic_flux_density(core, magnetic.get_coil(), &operatingPoint);
    auto current = operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current().value().get_waveform().value().get_data();
    auto [currentMinimum, currentMaximum] = extremes_abt1942(current);
    auto [fluxDensityMinimum, fluxDensityMaximum] = extremes_abt1942(fluxDensity.get_waveform().value().get_data());

    double gapReluctance = gap_reluctance_abt1942(core, InitialPermeability::get_initial_permeability(material, temperature, std::nullopt, frequency));
    double curvePeak = InitialPermeability::get_dc_operating_point_for_magnetomotive_force(
        material, numberTurns * currentMaximum, effectiveLength, effectiveArea, gapReluctance, temperature, frequency).second;
    INFO("i " << currentMinimum << " -> " << currentMaximum << " A; B " << fluxDensityMinimum << " -> " << fluxDensityMaximum
         << " T; curve peak " << curvePeak << " T; Lm " << resolve_dimensional_values(output.get_magnetizing_inductance()) * 1e6 << " uH");
    REQUIRE(currentMinimum == 0);
    REQUIRE_THAT(currentMaximum, Catch::Matchers::WithinRel(5.73, 0.01));
    CHECK_THAT(curvePeak, Catch::Matchers::WithinAbs(0.439, 0.005));
    CHECK_THAT(fluxDensityMinimum, Catch::Matchers::WithinAbs(0, 1e-9));
    CHECK_THAT(fluxDensityMaximum, Catch::Matchers::WithinRel(curvePeak, 1e-3));
    settings.reset();
}

TEST_CASE("Test_Dcm_Flyback_Flux_Density_Is_Faraday_From_Zero_With_The_Current_From_The_Curve", "[physical-model][magnetizing-inductance][dc-bias][abt-1942]") {
    // ABT #1942: slow_simulation's second operating point (flyback, material 95, 62:7:3, primary voltage given) has
    // its magnetizing current from the voltage: the flux swing is Faraday's, 0.464 T. The load line of that current's
    // average MMF lifted it onto a 0.588 T floor (0.588 -> 1.052 T, above saturation). The core comes back to zero MMF
    // in the dead time, so B is Faraday's anchored at zero there, and the magnetizing current is read back from the
    // material's curve on the gaps' load line: at the peak B_LL(N i_pk) is Faraday's peak.
    settings.reset();
    auto masJson = load_abt1942("slow_simulation.json");
    OpenMagnetics::Magnetic magnetic(masJson["magnetic"]);
    OpenMagnetics::Inputs inputs(masJson["inputs"]);
    REQUIRE(inputs.get_operating_points().size() > 1);
    auto operatingPoint = inputs.get_operating_points()[1];
    auto core = magnetic.get_core();
    auto material = core.resolve_material();
    double numberTurns = magnetic.get_coil().get_functional_description()[0].get_number_turns();
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    double effectiveLength = core.get_processed_description()->get_effective_parameters().get_effective_length();
    double temperature = operatingPoint.get_conditions().get_ambient_temperature();
    double frequency = operatingPoint.get_excitations_per_winding()[0].get_frequency();

    // Faraday on the primary voltage: the swing of the integral of v over N A_e.
    auto voltage = OpenMagnetics::Inputs::calculate_sampled_waveform(operatingPoint.get_excitations_per_winding()[0].get_voltage()->get_waveform().value(), frequency).get_data();
    double voltageAverage = 0;
    for (auto value : voltage) {
        voltageAverage += value / voltage.size();
    }
    double fluxLinkage = 0, fluxLinkageMinimum = 0, fluxLinkageMaximum = 0;
    for (auto value : voltage) {
        fluxLinkage += (value - voltageAverage) / frequency / voltage.size();
        fluxLinkageMinimum = std::min(fluxLinkageMinimum, fluxLinkage);
        fluxLinkageMaximum = std::max(fluxLinkageMaximum, fluxLinkage);
    }
    double faradaySwing = (fluxLinkageMaximum - fluxLinkageMinimum) / (numberTurns * effectiveArea);

    MagnetizingInductance magnetizingInductance;
    auto [output, fluxDensity] = magnetizingInductance.calculate_inductance_and_magnetic_flux_density(core, magnetic.get_coil(), &operatingPoint);
    auto [fluxDensityMinimum, fluxDensityMaximum] = extremes_abt1942(fluxDensity.get_waveform().value().get_data());
    auto current = operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current().value().get_waveform().value().get_data();
    auto [currentMinimum, currentMaximum] = extremes_abt1942(current);
    double gapReluctance = gap_reluctance_abt1942(core, InitialPermeability::get_initial_permeability(material, temperature, std::nullopt, frequency));
    double curveAtPeakCurrent = InitialPermeability::get_dc_operating_point_for_magnetomotive_force(
        material, numberTurns * currentMaximum, effectiveLength, effectiveArea, gapReluctance, temperature, frequency).second;
    INFO("B " << fluxDensityMinimum << " -> " << fluxDensityMaximum << " T, Faraday swing " << faradaySwing << " T; i_m "
         << currentMinimum << " -> " << currentMaximum << " A, B_LL(N i_pk) " << curveAtPeakCurrent << " T; Lm "
         << resolve_dimensional_values(output.get_magnetizing_inductance()) * 1e6 << " uH");
    REQUIRE_THAT(faradaySwing, Catch::Matchers::WithinAbs(0.464, 0.005));
    CHECK_THAT(fluxDensityMaximum - fluxDensityMinimum, Catch::Matchers::WithinRel(faradaySwing, 0.01));
    CHECK_THAT(fluxDensityMinimum, Catch::Matchers::WithinAbs(0, 0.002 * faradaySwing));
    CHECK_THAT(curveAtPeakCurrent, Catch::Matchers::WithinRel(fluxDensityMaximum, 1e-3));
    settings.reset();
}

namespace {
// low_filling_inductor's current lifted by `lift` amperes, through the magnetizing model: the small-signal model's
// signature (the average of B is the load line's B at the average MMF, the ripple linear at the incremental
// permeability of that point, whose inductance the model reports), and the curve-from-zero model's swing, which must
// differ from it for the case to tell the two apart. `reverses` says which side of zero the lifted current's minimum
// is required on: above it (CCM), or below it (a current reversing through zero; an imposed current reverses by any
// amount above rounding, ABT #1942).
void check_small_signal_load_line_abt1942(double lift, bool reverses) {
    settings.reset();
    auto masJson = load_abt1942("low_filling_inductor.json");
    auto& currentData = masJson["inputs"]["operatingPoints"][0]["excitationsPerWinding"][0]["current"]["waveform"]["data"];
    for (auto& datum : currentData) {
        datum = datum.get<double>() + lift;
    }
    OpenMagnetics::Magnetic magnetic(masJson["magnetic"]);
    OpenMagnetics::Inputs inputs(masJson["inputs"]);
    auto operatingPoint = inputs.get_operating_points()[0];
    auto core = magnetic.get_core();
    auto material = core.resolve_material();
    double numberTurns = magnetic.get_coil().get_functional_description()[0].get_number_turns();
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    double effectiveLength = core.get_processed_description()->get_effective_parameters().get_effective_length();
    double temperature = operatingPoint.get_conditions().get_ambient_temperature();
    double frequency = operatingPoint.get_excitations_per_winding()[0].get_frequency();

    MagnetizingInductance magnetizingInductance;
    auto [output, fluxDensity] = magnetizingInductance.calculate_inductance_and_magnetic_flux_density(core, magnetic.get_coil(), &operatingPoint);
    double inductance = resolve_dimensional_values(output.get_magnetizing_inductance());
    auto currentWaveform = operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current().value().get_waveform().value();
    auto [currentMinimum, currentMaximum] = extremes_abt1942(currentWaveform.get_data());
    double currentAverage = OpenMagnetics::Inputs::calculate_waveform_average(currentWaveform);
    auto fluxDensityWaveform = fluxDensity.get_waveform().value();
    auto [fluxDensityMinimum, fluxDensityMaximum] = extremes_abt1942(fluxDensityWaveform.get_data());
    double fluxDensityAverage = OpenMagnetics::Inputs::calculate_waveform_average(fluxDensityWaveform);

    double gapReluctance = gap_reluctance_abt1942(core, InitialPermeability::get_initial_permeability(material, temperature, std::nullopt, frequency));
    auto curve = [&](double current) {
        double magnetomotiveForce = numberTurns * current;
        return std::copysign(InitialPermeability::get_dc_operating_point_for_magnetomotive_force(
            material, std::fabs(magnetomotiveForce), effectiveLength, effectiveArea, gapReluctance, temperature, frequency).second, magnetomotiveForce);
    };
    double loadLineFluxDensity = curve(currentAverage);
    double incrementalSwing = inductance * (currentMaximum - currentMinimum) / (numberTurns * effectiveArea);
    double curveSwing = curve(currentMaximum) - curve(currentMinimum);
    INFO("lift " << lift << " A: i " << currentMinimum << " -> " << currentMaximum << " A (avg " << currentAverage << "); B " << fluxDensityMinimum << " -> " << fluxDensityMaximum << " T (avg " << fluxDensityAverage
         << "); load line B(N I_avg) " << loadLineFluxDensity << " T; swing incremental " << incrementalSwing << " T, curve "
         << curveSwing << " T; Lm " << inductance * 1e6 << " uH");
    if (reverses) {
        REQUIRE(currentMinimum < -0.1 * currentMaximum);
    }
    else {
        REQUIRE(currentMinimum > 0);
    }
    // The case tells the two models apart only if their swings differ.
    REQUIRE(std::fabs(curveSwing - incrementalSwing) > 0.05 * incrementalSwing);
    CHECK_THAT(fluxDensityMaximum - fluxDensityMinimum, Catch::Matchers::WithinRel(incrementalSwing, 2e-3));
    if (!reverses) {
        // A DC point the current rides on: the offset is the load line's. Through zero the waveform average is not
        // what today's path anchors (its harmonic DC), which ABT #1959 revisits; only its swing is pinned here.
        CHECK_THAT(fluxDensityAverage, Catch::Matchers::WithinRel(loadLineFluxDensity, 2e-3));
    }
    settings.reset();
}
}

TEST_CASE("Test_Ccm_Choke_Keeps_The_Small_Signal_Load_Line_Flux_Density", "[physical-model][magnetizing-inductance][dc-bias][abt-1942]") {
    // ABT #1942: the curve-from-zero path is for a magnetizing current that comes back to zero. A CCM choke (here
    // low_filling_inductor's current lifted by 3 A, so it never reaches zero) rides a minor loop on its DC point and
    // keeps the small-signal model.
    check_small_signal_load_line_abt1942(3.0, false);
}

TEST_CASE("Test_Magnetizing_Current_Through_Zero_Keeps_The_Small_Signal_Load_Line_Flux_Density", "[physical-model][magnetizing-inductance][dc-bias][abt-1942]") {
    // ABT #1942, unipolar DCM only (Alf, 2026-10-18): a magnetizing current that goes on through zero to the other
    // polarity (here low_filling_inductor's 0 -> 5.73 A lowered by 2 A, so it swings -2 -> 3.73 A) is not a core
    // starting from zero every cycle and keeps today's small-signal model (its own model is ABT #1959). This current
    // falls by its whole 5.73 A in one sample, so a tolerance of one sample step would have let it through: an imposed
    // current reaches zero only within rounding.
    check_small_signal_load_line_abt1942(-2.0, true);
}

TEST_CASE("Test_Magnetomotive_Force_For_Flux_Density_Throws_Past_The_Curve", "[physical-model][permeability][dc-bias][abt-1942]") {
    // ABT #1942 (Alf, 2026-10-18, no silent extrapolation): the inverse load line reads the magnetisation curve and
    // nothing beyond it. On the table it inverts the forward crossing; past the table's last point (1000 Hsat) it
    // throws instead of going on along the vacuum slope: a SaturationException above the material's B_sat, and a
    // data error when the curve ends below B_sat.
    settings.reset();
    auto material = Core::resolve_material(std::string("3C95"));
    auto curve = InitialPermeability::get_magnetisation_curve(material, 25, std::nullopt);
    const double effectiveLength = 0.05, effectiveArea = 1e-4, gapReluctance = 1e6;
    double lastFluxDensity = curve->fluxDensity.back();
    INFO("3C95 at 25 C: curve ends at " << lastFluxDensity << " T, H " << curve->fieldStrength.back() << " A/m; B_sat " << curve->saturationFluxDensity << " T");
    REQUIRE(lastFluxDensity > curve->saturationFluxDensity);
    for (double fraction : {0.3, 0.999}) {
        auto [fieldStrength, magnetomotiveForce] = InitialPermeability::get_magnetomotive_force_for_flux_density(
            *curve, fraction * lastFluxDensity, effectiveLength, effectiveArea, gapReluctance);
        auto forward = InitialPermeability::get_dc_operating_point_for_magnetomotive_force(*curve, magnetomotiveForce, effectiveLength, effectiveArea, gapReluctance);
        CHECK_THAT(forward.second, Catch::Matchers::WithinRel(fraction * lastFluxDensity, 1e-9));
        CHECK_THAT(forward.first, Catch::Matchers::WithinRel(fieldStrength, 1e-9));
    }
    CHECK_THROWS_AS(InitialPermeability::get_magnetomotive_force_for_flux_density(
        *curve, 1.001 * lastFluxDensity, effectiveLength, effectiveArea, gapReluctance), SaturationException);

    // A curve that ends below its material's B_sat: past its end is a contradiction in the data, not saturation.
    InitialPermeability::MagnetisationCurve shortCurve{{0, 10, 20}, {0, 0.1, 0.2}, "short", 0.5, false};
    bool saturationThrown = false, otherThrown = false;
    try {
        InitialPermeability::get_magnetomotive_force_for_flux_density(shortCurve, 0.3, effectiveLength, effectiveArea, gapReluctance);
    }
    catch (const SaturationException&) {
        saturationThrown = true;
    }
    catch (const std::exception&) {
        otherThrown = true;
    }
    CHECK_FALSE(saturationThrown);
    CHECK(otherThrown);
    settings.reset();
}

TEST_CASE("Test_Dcm_Imposed_Current_Past_The_Knee_Is_Flagged_On_Hsat", "[physical-model][magnetizing-inductance][dc-bias][saturation][abt-1942]") {
    // ABT #1942 (Alf, 2026-10-18, peak H against the material's Hsat): where B is read on the magnetisation curve it
    // rises past the knee only on the vacuum slope, so it can stay under B_sat while the core is driven deep into
    // saturation. low_filling_inductor (Fair-Rite 78, 0 -> 5.73 A) peaks at 0.439 T on its curve, under B_sat 0.455 T,
    // but the curve gives it that flux at a field strength above Hsat 397 A/m (B(Hsat) = 0.404 T): the saturation check
    // analyze_results and the advisers make flags it, on the field strength (813 A/m). At 40 % of the current it stays
    // under Hsat and is not flagged.
    settings.reset();
    auto verdict_at = [](double scale) {
        auto masJson = load_abt1942("low_filling_inductor.json");
        auto& currentData = masJson["inputs"]["operatingPoints"][0]["excitationsPerWinding"][0]["current"]["waveform"]["data"];
        for (auto& datum : currentData) {
            datum = datum.get<double>() * scale;
        }
        OpenMagnetics::Magnetic magnetic(masJson["magnetic"]);
        OpenMagnetics::Inputs inputs(masJson["inputs"]);
        auto operatingPoint = inputs.get_operating_points()[0];
        auto core = magnetic.get_core();
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();
        MagnetizingInductance magnetizingInductance;
        auto fluxDensity = magnetizingInductance.calculate_inductance_and_magnetic_flux_density(core, magnetic.get_coil(), &operatingPoint).second;
        auto dcBiasPoint = magnetizingInductance.get_last_dc_bias_point();
        REQUIRE(dcBiasPoint);
        REQUIRE(dcBiasPoint->peakMagneticFieldStrength);
        return ResultsAnalyzer::check_saturation(core.resolve_material(), ResultsAnalyzer::get_peak_flux_density(fluxDensity), temperature, 1.0, dcBiasPoint);
    };

    auto full = verdict_at(1.0);
    INFO(full.message);
    REQUIRE(full.peakFieldStrength);
    REQUIRE(full.saturationFieldStrength);
    CHECK(full.peakFieldStrength.value() > full.saturationFieldStrength.value());
    CHECK_FALSE(full.fluxDensityAboveSaturation);
    CHECK(full.fieldStrengthAboveSaturation);
    CHECK(full.saturates);

    auto reduced = verdict_at(0.4);
    INFO(reduced.message);
    REQUIRE(reduced.peakFieldStrength);
    CHECK(reduced.peakFieldStrength.value() < reduced.saturationFieldStrength.value());
    CHECK_FALSE(reduced.saturates);
    settings.reset();
}

TEST_CASE("Test_Dcm_Voltage_Driven_Paramagnetic_Core_Reads_The_Current_On_Mu0", "[physical-model][magnetizing-inductance][dc-bias][abt-1942]") {
    // ABT #1942 (Alf, 2026-10-19): a core the loss solve holds at its Curie temperature is magnetized as the
    // paramagnetic core it is above it (mu_r = 1, MagneticSimulator::make_paramagnetic_magnetic). Its magnetisation
    // curve is mu0*H at every H; the table that stores it ends at 1000 Hsat(T), a few A/m and a few mT just under the
    // Curie point (Hsat(T) goes to zero there). A voltage-driven DCM point reads its magnetizing current back from the
    // curve: on the paramagnetic core that is mu0's, H = B/mu0, at any B. Before, the readback refused a Faraday flux
    // past the table's end ("integrates to only 0.006346 T at 5049.8 A/m"), and the MagneticAdviser lost the whole call
    // on slow_simulation, error_paint_field and catalogue_missing_turns. slow_simulation's second operating point
    // (flyback, material 95, primary voltage given), held 1 K under 95's Curie temperature.
    settings.reset();
    auto masJson = load_abt1942("slow_simulation.json");
    OpenMagnetics::Magnetic magnetic(masJson["magnetic"]);
    OpenMagnetics::Inputs inputs(masJson["inputs"]);
    REQUIRE(inputs.get_operating_points().size() > 1);
    auto operatingPoint = inputs.get_operating_points()[1];
    auto paramagnetic = MagneticSimulator::make_paramagnetic_magnetic(magnetic);
    auto core = paramagnetic.get_core();
    auto material = core.resolve_material();
    auto curieTemperature = magnetic.get_core().resolve_material().get_curie_temperature();
    REQUIRE(curieTemperature);
    double temperature = curieTemperature.value() - 1;
    operatingPoint.get_mutable_conditions().set_ambient_temperature(temperature);
    double numberTurns = magnetic.get_coil().get_functional_description()[0].get_number_turns();
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    double frequency = operatingPoint.get_excitations_per_winding()[0].get_frequency();

    REQUIRE(InitialPermeability::is_paramagnetic(material, temperature));
    REQUIRE_FALSE(InitialPermeability::is_paramagnetic(magnetic.get_core().resolve_material(), temperature));
    auto curve = InitialPermeability::get_magnetisation_curve(material, temperature, frequency);
    REQUIRE(curve->paramagnetic);

    MagnetizingInductance magnetizingInductance;
    auto [output, fluxDensity] = magnetizingInductance.calculate_inductance_and_magnetic_flux_density(core, paramagnetic.get_coil(), &operatingPoint);
    auto [fluxDensityMinimum, fluxDensityMaximum] = extremes_abt1942(fluxDensity.get_waveform().value().get_data());
    auto current = operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current().value().get_waveform().value().get_data();
    auto [currentMinimum, currentMaximum] = extremes_abt1942(current);
    // mu_r = 1 in the core: N i = B A_e R_total(mu_r = 1) at every sample, the peak included.
    double totalReluctance = ReluctanceModel::factory()->get_core_reluctance(core, 1.0).get_core_reluctance();
    double expectedPeakCurrent = fluxDensityMaximum * effectiveArea * totalReluctance / numberTurns;
    INFO("paramagnetic 95 at " << temperature << " C: curve table ends at " << curve->fluxDensity.back() << " T, "
         << curve->fieldStrength.back() << " A/m; B " << fluxDensityMinimum << " -> " << fluxDensityMaximum << " T; i_m "
         << currentMinimum << " -> " << currentMaximum << " A, mu0 load line " << expectedPeakCurrent << " A");
    REQUIRE(fluxDensityMaximum > curve->fluxDensity.back());
    CHECK_THAT(fluxDensityMinimum, Catch::Matchers::WithinAbs(0, 0.002 * fluxDensityMaximum));
    CHECK_THAT(currentMaximum, Catch::Matchers::WithinRel(expectedPeakCurrent, 1e-3));

    // The inverse itself, far past the table: H = B/mu0.
    auto [fieldStrength, magnetomotiveForce] = InitialPermeability::get_magnetomotive_force_for_flux_density(*curve, 0.3, 0.05, 1e-4, 1e6);
    CHECK_THAT(fieldStrength, Catch::Matchers::WithinRel(0.3 / constants.vacuumPermeability, 1e-12));
    CHECK_THAT(magnetomotiveForce, Catch::Matchers::WithinRel(fieldStrength * 0.05 + 0.3 * 1e-4 * 1e6, 1e-12));
    settings.reset();
}

namespace {
// The web DMC wizard's switching-ripple record (ABT #1764/#1786): two windings, each with the 50 Hz line current
// (4.8 A peak) and a 0.44 A, 100 kHz ripple, 16384 samples; the currents scaled by currentScale.
std::pair<OpenMagnetics::Inputs, OperatingPoint> load_dmc_line_cycle_abt1786(double currentScale) {
    std::ifstream file(OpenMagneticsTesting::get_test_data_path(std::source_location::current(),
                                                                "web_adviser_flows/dmc_calculate_advised_magnetics_switching_ripple.json"));
    REQUIRE(file.good());
    json inputsJson = json::parse(file).at("inputs");
    for (auto& excitation : inputsJson["operatingPoints"][0]["excitationsPerWinding"]) {
        excitation.erase("magnetizingCurrent");
        for (auto& datum : excitation["current"]["waveform"]["data"]) {
            datum = datum.get<double>() * currentScale;
        }
    }
    OpenMagnetics::Inputs inputs(inputsJson);
    return {inputs, inputs.get_operating_point(0)};
}
}

TEST_CASE("Test_Line_Cycle_Swing_Follows_The_Load_Line_On_The_Magnetisation_Curve", "[physical-model][magnetizing-inductance][dc-bias][line-cycle][abt-1786]") {
    // ABT #1786: a line-cycle record under an imposed MMF (a DMC choke carrying its 50 Hz line current with the
    // switching ripple) read each switching cycle's line flux density linearly, N i / (R(mu_peak) A_e), with mu_peak the
    // reversible permeability at the line peak: 0.354 T at the peak on the web DMC fixture's ungapped T 48/23/19 MPP 300
    // with 30 + 30 turns (H_pk 2381 A/m). The material's magnetisation curve, B(H) = mu0 * integral of mu_rev, lies well
    // above mu0 mu_rev(H) H there. Alf, 2026-10-18: the line swing follows the load line on the curve like the DC case,
    // so the flux density at the line peak is B_LL(N i_pk) (each cycle at its own mean; the ripple about it is small),
    // the most negative line cycle its mirror, and the DC point carries the load line's H at the peak for the H_sat check.
    settings.reset();
    auto [inputs, operatingPoint] = load_dmc_line_cycle_abt1786(1.0);
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic("T 48/23/19", json::array(), {30, 30}, 1, "MPP 300");
    auto core = magnetic.get_core();
    auto material = core.resolve_material();
    double numberTurns = 30;
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    double effectiveLength = core.get_processed_description()->get_effective_parameters().get_effective_length();
    double temperature = operatingPoint.get_conditions().get_ambient_temperature();
    double frequency = operatingPoint.get_excitations_per_winding()[0].get_frequency();

    MagnetizingInductance magnetizingInductance;
    magnetizingInductance.set_topology(inputs.get_design_requirements().get_topology());
    auto [output, fluxDensity] = magnetizingInductance.calculate_inductance_and_magnetic_flux_density(core, magnetic.get_coil(), &operatingPoint);
    auto current = operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current().value().get_waveform().value().get_data();
    auto [currentMinimum, currentMaximum] = extremes_abt1942(current);
    auto [fluxDensityMinimum, fluxDensityMaximum] = extremes_abt1942(fluxDensity.get_waveform().value().get_data());

    double gapReluctance = gap_reluctance_abt1942(core, InitialPermeability::get_initial_permeability(material, temperature, std::nullopt, frequency));
    auto [peakFieldStrength, curvePeak] = InitialPermeability::get_dc_operating_point_for_magnetomotive_force(
        material, numberTurns * currentMaximum, effectiveLength, effectiveArea, gapReluctance, temperature, frequency);
    double curveValley = InitialPermeability::get_dc_operating_point_for_magnetomotive_force(
        material, numberTurns * fabs(currentMinimum), effectiveLength, effectiveArea, gapReluctance, temperature, frequency).second;
    auto dcBiasPoint = magnetizingInductance.get_last_dc_bias_point();
    INFO("i_m " << currentMinimum << " -> " << currentMaximum << " A; B " << fluxDensityMinimum << " -> " << fluxDensityMaximum
         << " T; curve " << -curveValley << " -> " << curvePeak << " T at H " << peakFieldStrength << " A/m");
    REQUIRE(currentMaximum > 9);
    REQUIRE(currentMinimum < -9);
    CHECK(curvePeak > 0.5);
    CHECK_THAT(fluxDensityMaximum, Catch::Matchers::WithinRel(curvePeak, 0.03));
    CHECK_THAT(fluxDensityMinimum, Catch::Matchers::WithinRel(-curveValley, 0.03));
    REQUIRE(dcBiasPoint);
    REQUIRE(dcBiasPoint->peakMagneticFieldStrength);
    // At the record's largest |N i|, which is the negative line peak here.
    double largestFieldStrength = InitialPermeability::get_dc_operating_point_for_magnetomotive_force(
        material, numberTurns * std::max(currentMaximum, fabs(currentMinimum)), effectiveLength, effectiveArea, gapReluctance, temperature, frequency).first;
    CHECK_THAT(dcBiasPoint->peakMagneticFieldStrength.value(), Catch::Matchers::WithinRel(largestFieldStrength, 1e-6));
    settings.reset();
}

TEST_CASE("Test_Line_Cycle_Swing_Past_The_Curve_End_Throws", "[physical-model][magnetizing-inductance][dc-bias][line-cycle][abt-1786]") {
    // ABT #1786: a line MMF past the load line at the end of the magnetisation curve (1000 times H_sat) has no flux
    // density on it; the line swing throws there rather than extrapolate (as the inverse does, ABT #1942). The DMC
    // fixture's currents scaled 10^5 drive the MPP 300 toroid's line peak to ~2.4e8 A/m.
    settings.reset();
    auto [inputs, operatingPoint] = load_dmc_line_cycle_abt1786(1e5);
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic("T 48/23/19", json::array(), {30, 30}, 1, "MPP 300");
    MagnetizingInductance magnetizingInductance;
    magnetizingInductance.set_topology(inputs.get_design_requirements().get_topology());
    CHECK_THROWS_AS(magnetizingInductance.calculate_inductance_and_magnetic_flux_density(magnetic.get_core(), magnetic.get_coil(), &operatingPoint),
                    SaturationException);
    settings.reset();
}


// ABT #2002 (Alf, 2026-10-24): the winding's resistive drop is taken out only to judge volt-second balance. The
// magnetizing current is derived from the voltage as given, which MKF treats as the induced voltage (converter
// models and the web give ideal voltages): a wound primary (R0 > 0) and the same coil unwound (no length, no
// resistance) give the same magnetizing current from the same voltage, even one carrying R0 i0(t). Taking the
// drop out here tilted a DCM flyback's dead time and gave a common-mode choke's 50 Hz line point a flux of -R i.
TEST_CASE("Test_Magnetizing_Current_Is_Derived_From_The_Voltage_As_Given", "[physical-model][magnetizing-inductance][abt-2002]") {
    settings.reset();
    auto path = std::filesystem::path{std::source_location::current().file_name()}
                    .parent_path().append("..").append("MAS").append("examples").append("02_flyback_efd25_3c95.json");
    REQUIRE(std::filesystem::exists(path));
    auto mas = OpenMagneticsTesting::mas_loader(path.string());
    auto wound = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    if (!wound.get_coil().get_turns_description()) {
        wound.get_mutable_coil().wind();
    }
    REQUIRE(wound.get_coil().get_turns_description());
    OpenMagnetics::Magnetic unwound = wound;
    unwound.get_mutable_coil().set_turns_description(std::nullopt);

    OperatingPoint given = mas.get_inputs().get_operating_points()[0];
    REQUIRE(given.get_excitations_per_winding().size() == 2);
    double temperature = given.get_conditions().get_ambient_temperature();
    double primaryResistance = OpenMagnetics::Inputs::get_winding_dc_resistances(wound.get_coil(), temperature).value()[0];
    REQUIRE(primaryResistance > 0);

    // Both operating points on one sampled grid, with the same primary voltage carrying R0 i0(t).
    auto primary = given.get_excitations_per_winding()[0];
    double frequency = primary.get_frequency();
    auto sampledVoltage = OpenMagnetics::Inputs::calculate_sampled_waveform(
        OpenMagnetics::Inputs::standardize_waveform(primary.get_voltage().value(), frequency).get_waveform().value(), frequency);
    auto sampledCurrent = OpenMagnetics::Inputs::calculate_sampled_waveform(
        OpenMagnetics::Inputs::standardize_waveform(primary.get_current().value(), frequency).get_waveform().value(), frequency);
    REQUIRE(sampledVoltage.get_data().size() == sampledCurrent.get_data().size());
    Waveform voltageWithDrop = sampledVoltage;
    std::vector<double> dataWithDrop;
    for (size_t index = 0; index < sampledVoltage.get_data().size(); ++index) {
        dataWithDrop.push_back(sampledVoltage.get_data()[index] + primaryResistance * sampledCurrent.get_data()[index]);
    }
    voltageWithDrop.set_data(dataWithDrop);

    auto with_primary = [&](const Waveform& voltage) {
        OperatingPoint operatingPoint = given;
        for (auto& excitation : operatingPoint.get_mutable_excitations_per_winding()) {
            excitation.set_magnetizing_current(std::nullopt);
            excitation.set_magnetic_flux_density(std::nullopt);
            excitation.set_magnetic_field_strength(std::nullopt);
        }
        SignalDescriptor voltageSignal;
        voltageSignal.set_waveform(voltage);
        SignalDescriptor currentSignal;
        currentSignal.set_waveform(sampledCurrent);
        operatingPoint.get_mutable_excitations_per_winding()[0].set_voltage(voltageSignal);
        operatingPoint.get_mutable_excitations_per_winding()[0].set_current(currentSignal);
        return operatingPoint;
    };
    OperatingPoint woundOperatingPoint = with_primary(voltageWithDrop);
    OperatingPoint unwoundOperatingPoint = with_primary(voltageWithDrop);

    MagnetizingInductance model;
    auto woundResult = model.calculate_inductance_and_magnetic_flux_density(wound, &woundOperatingPoint);
    auto unwoundResult = model.calculate_inductance_and_magnetic_flux_density(unwound, &unwoundOperatingPoint);
    REQUIRE(resolve_dimensional_values(woundResult.first.get_magnetizing_inductance()) ==
            resolve_dimensional_values(unwoundResult.first.get_magnetizing_inductance()));

    auto woundCurrent = woundOperatingPoint.get_excitations_per_winding()[0].get_magnetizing_current()->get_waveform()->get_data();
    auto unwoundCurrent = unwoundOperatingPoint.get_excitations_per_winding()[0].get_magnetizing_current()->get_waveform()->get_data();
    REQUIRE(woundCurrent.size() == unwoundCurrent.size());
    double peakToPeak = *std::max_element(unwoundCurrent.begin(), unwoundCurrent.end()) - *std::min_element(unwoundCurrent.begin(), unwoundCurrent.end());
    REQUIRE(peakToPeak > 0);
    double largestDifference = 0;
    for (size_t index = 0; index < woundCurrent.size(); ++index) {
        largestDifference = std::max(largestDifference, std::fabs(woundCurrent[index] - unwoundCurrent[index]));
    }
    INFO("R0 = " << primaryResistance << " Ohm, magnetizing current pk-pk " << peakToPeak << " A, largest difference " << largestDifference << " A");
    CHECK(largestDifference < 1e-9 * peakToPeak);
}
