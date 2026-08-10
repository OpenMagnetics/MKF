#include <source_location>
#include "processors/CircuitSimulatorInterface.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/Reluctance.h"
#include "constructive_models/Bobbin.h"
#include "constructive_models/Magnetic.h"
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
            "open", "name": "ETD 54/28/19", "type": "standard"}, "type": "twoPieceSet"},
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
        return json{{"type", "magneticEpoxy"}, {"thickness", 0.0001}, {"material", materialName}};
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
        buildCore(json{{"type", "epoxy"}, {"thickness", 0.0001}}), numberTurns, 25));  // non-magnetic coating
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
