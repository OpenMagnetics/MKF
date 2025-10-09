#include "advisers/CoreCrossReferencer.h"
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

using namespace MAS;
using namespace OpenMagnetics;


SUITE(CoreCrossReferencer) {

    TEST(Test_All_Core_Materials) {
        settings->reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreCrossReferencer coreCrossReferencer;

        std::string coreName = "EC 35/17/10 - 3C91 - Gapped 1.000 mm";
        Core core = find_core_by_name(coreName);

        double temperature = 20;
        auto label = WaveformLabel::TRIANGULAR;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double frequency = 100000;
        double magnetizingInductance = 100e-6;
        int64_t numberTurns = 28;

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              label,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              offset);


        auto crossReferencedCores = coreCrossReferencer.get_cross_referenced_core(core, numberTurns, inputs, 5);


        CHECK(crossReferencedCores.size() > 0);
        CHECK(crossReferencedCores[0].first.get_name() == "EC 35/17/10 - 3C94 - Gapped 1.000 mm");
    }

    TEST(Test_All_Core_Materials_Same_Material) {
        settings->reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreCrossReferencer coreCrossReferencer;
        coreCrossReferencer.use_only_reference_material(true);

        std::string coreName = "EC 35/17/10 - 3C91 - Gapped 1.000 mm";
        Core core = find_core_by_name(coreName);

        double temperature = 20;
        auto label = WaveformLabel::TRIANGULAR;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double frequency = 100000;
        double magnetizingInductance = 100e-6;
        int64_t numberTurns = 28;

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              label,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              offset);


        auto crossReferencedCores = coreCrossReferencer.get_cross_referenced_core(core, numberTurns, inputs, 5);


        CHECK(crossReferencedCores.size() > 0);
        CHECK(crossReferencedCores[0].first.get_name() == "EP 20 - 3C91 - Gapped 0.605 mm");
    }

    TEST(Test_All_Core_Materials_Same_Material_Maximum_Height) {
        settings->reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreCrossReferencer coreCrossReferencer;
        coreCrossReferencer.use_only_reference_material(true);

        std::string coreName = "EC 35/17/10 - 3C91 - Gapped 1.000 mm";
        Core core = find_core_by_name(coreName);

        double temperature = 20;
        auto label = WaveformLabel::TRIANGULAR;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double frequency = 100000;
        double magnetizingInductance = 100e-6;
        int64_t numberTurns = 28;

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              label,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              offset);

        MaximumDimensions maximumDimensions;
        maximumDimensions.set_depth(0.04);
        maximumDimensions.set_height(0.025);
        maximumDimensions.set_width(0.025);
        inputs.get_mutable_design_requirements().set_maximum_dimensions(maximumDimensions);


        auto crossReferencedCores = coreCrossReferencer.get_cross_referenced_core(core, numberTurns, inputs, 5);


        CHECK(crossReferencedCores.size() > 0);
        CHECK(crossReferencedCores[0].first.get_name() == "EP 20 - 3C91 - Gapped 0.605 mm");
    }

    TEST(Test_All_Core_Materials_Only_TDK) {
        settings->reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreCrossReferencer coreCrossReferencer;
        coreCrossReferencer.use_only_manufacturer("TDK");

        std::string coreName = "EC 35/17/10 - 3C91 - Gapped 1.000 mm";
        Core core = find_core_by_name(coreName);

        double temperature = 20;
        auto label = WaveformLabel::TRIANGULAR;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double frequency = 100000;
        double magnetizingInductance = 100e-6   ;
        int64_t numberTurns = 28;

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                         magnetizingInductance,
                                                                                         temperature,
                                                                                         label,
                                                                                         peakToPeak,
                                                                                         dutyCycle,
                                                                                         offset);


        auto crossReferencedCores = coreCrossReferencer.get_cross_referenced_core(core, numberTurns, inputs, 5);


        CHECK(crossReferencedCores.size() > 0);
        CHECK(crossReferencedCores[0].first.get_name() == "ETD 34/17/11 - N87 - Gapped 2.500 mm");
    }

    TEST(Test_All_Core_Materials_Powder) {
        settings->reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreCrossReferencer coreCrossReferencer;

        std::string coreName = "E 25/9.5/6.3 - XFlux 60 - Ungapped";
        Core core = find_core_by_name(coreName);

        double temperature = 20;
        auto label = WaveformLabel::TRIANGULAR;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double frequency = 100000;
        double magnetizingInductance = 100e-6;  
        int64_t numberTurns = 28;

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                         magnetizingInductance,
                                                                                         temperature,
                                                                                         label,
                                                                                         peakToPeak,
                                                                                         dutyCycle,
                                                                                         offset);


        auto crossReferencedCores = coreCrossReferencer.get_cross_referenced_core(core, numberTurns, inputs, 5);


        CHECK(crossReferencedCores.size() > 0);
        CHECK(crossReferencedCores[0].first.get_name() == "E 25/9.5/6.3 - Kool Mµ Hƒ 60 - Ungapped");
    }

    TEST(Test_All_Core_Materials_Only_Micrometals) {
        settings->reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreCrossReferencer coreCrossReferencer;
        coreCrossReferencer.use_only_manufacturer("Micrometals");

        std::string coreName = "E 25/9.5/6.3 - XFlux 60 - Ungapped";
        Core core = find_core_by_name(coreName);

        double temperature = 20;
        auto label = WaveformLabel::TRIANGULAR;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double frequency = 100000;
        double magnetizingInductance = 100e-6;
        int64_t numberTurns = 28;

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                         magnetizingInductance,
                                                                                         temperature,
                                                                                         label,
                                                                                         peakToPeak,
                                                                                         dutyCycle,
                                                                                         offset);


        auto crossReferencedCores = coreCrossReferencer.get_cross_referenced_core(core, numberTurns, inputs, 5);


        CHECK(crossReferencedCores.size() > 0);
        CHECK(crossReferencedCores[0].first.get_name() == "E 25/9.5/6.3 - FS 26 - Ungapped");
    }

    TEST(Test_Cross_Reference_Core_Web_0) {
        settings->reset();
        settings->set_use_only_cores_in_stock(false);
        clear_databases();
        OperatingPoint operatingPoint;
        CoreCrossReferencer coreCrossReferencer;

        std::string shapeName = "PQ 40/40";
        std::string materialName = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, materialName);


        double temperature = 20;
        auto label = WaveformLabel::TRIANGULAR;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double frequency = 100000;
        double magnetizingInductance = 100e-6;
        int64_t numberTurns = 16;

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                         magnetizingInductance,
                                                                                         temperature,
                                                                                         label,
                                                                                         peakToPeak,
                                                                                         dutyCycle,
                                                                                         offset);


        auto crossReferencedCores = coreCrossReferencer.get_cross_referenced_core(core, numberTurns, inputs, 20);

        CHECK(crossReferencedCores.size() > 0);
        CHECK(crossReferencedCores[0].first.get_name() == "EQ 42/20/28 - Kool Mµ MAX 26 - Ungapped");

        auto scorings = coreCrossReferencer.get_scorings();
        auto scoredValues = coreCrossReferencer.get_scored_values();
        json results;
        results["cores"] = json::array();
        results["scorings"] = json::array();
        for (auto& [core, scoring] : crossReferencedCores) {
            std::string name = core.get_name().value();

            json coreJson;
            to_json(coreJson, core);
            results["cores"].push_back(coreJson);
            results["scorings"].push_back(scoring);

            json result;
            result["scoringPerFilter"] = json();
            result["scoredValuePerFilter"] = json();
            for (auto& filter : magic_enum::enum_names<CoreCrossReferencerFilters>()) {
                std::string filterString(filter);

                result["scoringPerFilter"][filterString] = scorings[name][magic_enum::enum_cast<CoreCrossReferencerFilters>(filterString).value()];
                result["scoredValuePerFilter"][filterString] = scoredValues[name][magic_enum::enum_cast<CoreCrossReferencerFilters>(filterString).value()];
                CHECK(!std::isnan(scorings[name][magic_enum::enum_cast<CoreCrossReferencerFilters>(filterString).value()]));
                CHECK(!std::isnan(scoredValues[name][magic_enum::enum_cast<CoreCrossReferencerFilters>(filterString).value()]));
            };
            results["data"].push_back(result);
        }

        settings->reset();
    }
}
