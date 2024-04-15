#include "Settings.h"
#include "CoreCrossReferencer.h"
#include "Utils.h"
#include "InputsWrapper.h"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>


SUITE(CoreCrossReferencer) {
    auto settings = OpenMagnetics::Settings::GetInstance();

    TEST(Test_All_Core_Materials) {
        settings->reset();
        OpenMagnetics::clear_databases();
        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreCrossReferencer coreCrossReferencer;

        std::string coreName = "EC 35/17/10 - 3C91 - Gapped 1.000 mm";
        OpenMagnetics::CoreWrapper core = OpenMagnetics::find_core_by_name(coreName);

        double temperature = 20;
        auto label = OpenMagnetics::WaveformLabel::TRIANGULAR;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double frequency = 100000;
        double magnetizingInductance = 100e-6   ;
        int64_t numberTurns = 28;

        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                         magnetizingInductance,
                                                                                         temperature,
                                                                                         label,
                                                                                         peakToPeak,
                                                                                         dutyCycle,
                                                                                         offset);


        auto crossReferencedCores = coreCrossReferencer.get_cross_referenced_core(core, numberTurns, inputs, 5);


        CHECK(crossReferencedCores.size() > 0);
        CHECK(crossReferencedCores[0].first.get_name() == "ETD 29/16/10 - 3C97 - Gapped 1.000 mm");
        // for (auto [core, scoring] : crossReferencedCores) {
            // std::cout << core.get_name().value() << std::endl;
        // }
    }

    TEST(Test_All_Core_Materials_Only_TDK) {
        settings->reset();
        OpenMagnetics::clear_databases();
        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreCrossReferencer coreCrossReferencer;
        coreCrossReferencer.use_only_manufacturer("TDK");

        std::string coreName = "EC 35/17/10 - 3C91 - Gapped 1.000 mm";
        OpenMagnetics::CoreWrapper core = OpenMagnetics::find_core_by_name(coreName);

        double temperature = 20;
        auto label = OpenMagnetics::WaveformLabel::TRIANGULAR;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double frequency = 100000;
        double magnetizingInductance = 100e-6   ;
        int64_t numberTurns = 28;

        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                         magnetizingInductance,
                                                                                         temperature,
                                                                                         label,
                                                                                         peakToPeak,
                                                                                         dutyCycle,
                                                                                         offset);


        auto crossReferencedCores = coreCrossReferencer.get_cross_referenced_core(core, numberTurns, inputs, 5);


        CHECK(crossReferencedCores.size() > 0);
        CHECK(crossReferencedCores[0].first.get_name() == "ER 42/22/15 - N27 - Gapped 1.000 mm");
    }

    TEST(Test_All_Core_Materials_Powder) {
        settings->reset();
        OpenMagnetics::clear_databases();
        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreCrossReferencer coreCrossReferencer;

        std::string coreName = "E 25/9.5/6.3 - XFlux 60 - Ungapped";
        OpenMagnetics::CoreWrapper core = OpenMagnetics::find_core_by_name(coreName);

        double temperature = 20;
        auto label = OpenMagnetics::WaveformLabel::TRIANGULAR;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double frequency = 100000;
        double magnetizingInductance = 100e-6;  
        int64_t numberTurns = 28;

        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                         magnetizingInductance,
                                                                                         temperature,
                                                                                         label,
                                                                                         peakToPeak,
                                                                                         dutyCycle,
                                                                                         offset);


        auto crossReferencedCores = coreCrossReferencer.get_cross_referenced_core(core, numberTurns, inputs, 5);


        CHECK(crossReferencedCores.size() > 0);
        CHECK(crossReferencedCores[0].first.get_name() == "E 25/9.5/6.3 - Kool Mµ 60 - Ungapped");
        for (auto [core, scoring] : crossReferencedCores) {
            std::cout << core.get_name().value() << std::endl;
        }
    }

    TEST(Test_All_Core_Materials_Only_Micrometals) {
        settings->reset();
        OpenMagnetics::clear_databases();
        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreCrossReferencer coreCrossReferencer;
        coreCrossReferencer.use_only_manufacturer("Micrometals");

        std::string coreName = "E 25/9.5/6.3 - XFlux 60 - Ungapped";
        OpenMagnetics::CoreWrapper core = OpenMagnetics::find_core_by_name(coreName);

        double temperature = 20;
        auto label = OpenMagnetics::WaveformLabel::TRIANGULAR;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double frequency = 100000;
        double magnetizingInductance = 100e-6;
        int64_t numberTurns = 28;

        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                         magnetizingInductance,
                                                                                         temperature,
                                                                                         label,
                                                                                         peakToPeak,
                                                                                         dutyCycle,
                                                                                         offset);


        auto crossReferencedCores = coreCrossReferencer.get_cross_referenced_core(core, numberTurns, inputs, 5);


        CHECK(crossReferencedCores.size() > 0);
        CHECK(crossReferencedCores[0].first.get_name() == "T 17.3/9.65/6.35 - parylene coated - OC 90 - Ungapped");
        for (auto [core, scoring] : crossReferencedCores) {
            std::cout << core.get_name().value() << std::endl;
        }
    }
}