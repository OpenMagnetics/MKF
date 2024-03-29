#include "CoreAdviser.h"
#include "WireAdviser.h"
#include "InputsWrapper.h"
#include "TestingUtils.h"

#include <UnitTest++.h>
#include <vector>

SUITE(WireAdviser) {
    OpenMagnetics::CoilFunctionalDescription coilFunctionalDescription;
    OpenMagnetics::Section section;
    OpenMagnetics::SignalDescriptor current;

    int64_t numberTurns = 42;
    double windingWindowWidth = 0.005;
    double windingWindowHeight = 0.015;
    double currentRms = 2;
    double currentEffectiveFrequency = 13456;
    double temperature = 22;
    uint8_t numberSections = 1;
    size_t maximumNumberResults = 1;


    void setup() {
        coilFunctionalDescription.set_isolation_side(OpenMagnetics::IsolationSide::PRIMARY);
        coilFunctionalDescription.set_name("primary");
        coilFunctionalDescription.set_number_parallels(1);
        coilFunctionalDescription.set_number_turns(numberTurns);
        coilFunctionalDescription.set_wire("Dummy");

        section.set_dimensions({windingWindowWidth, windingWindowHeight});

        OpenMagnetics::Processed processed;
        processed.set_rms(currentRms);
        processed.set_effective_frequency(currentEffectiveFrequency);
        current.set_processed(processed);
    }

    TEST(Test_Round) {
        numberTurns = 2;
        currentRms = 10;
        currentEffectiveFrequency = 134567;
        setup();
        OpenMagnetics::WireAdviser wireAdviser;
        wireAdviser.enableFoil(false);
        wireAdviser.enableRectangular(false);
        wireAdviser.enableLitz(false);
        wireAdviser.enableRound(true);
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription, section, current, temperature, numberSections, maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(OpenMagnetics::WireType::ROUND == OpenMagnetics::CoilWrapper::resolve_wire(masMagneticWithCoil).get_type());
    }

    TEST(Test_Litz) {
        numberTurns = 2;
        currentRms = 10;
        currentEffectiveFrequency = 134567;
        setup();
        OpenMagnetics::WireAdviser wireAdviser;
        wireAdviser.enableFoil(false);
        wireAdviser.enableRectangular(false);
        wireAdviser.enableLitz(true);
        wireAdviser.enableRound(false);
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription, section, current, temperature, numberSections, maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(OpenMagnetics::WireType::LITZ == OpenMagnetics::CoilWrapper::resolve_wire(masMagneticWithCoil).get_type());
    }

    TEST(Test_Rectangular) {
        numberTurns = 2;
        currentRms = 10;
        currentEffectiveFrequency = 134567;
        setup();
        OpenMagnetics::WireAdviser wireAdviser;
        wireAdviser.enableFoil(false);
        wireAdviser.enableRectangular(true);
        wireAdviser.enableLitz(false);
        wireAdviser.enableRound(false);
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription, section, current, temperature, numberSections, maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(OpenMagnetics::WireType::RECTANGULAR == OpenMagnetics::CoilWrapper::resolve_wire(masMagneticWithCoil).get_type());
    }

    TEST(Test_Foil) {
        numberTurns = 2;
        currentRms = 10;
        currentEffectiveFrequency = 134567;
        setup();
        OpenMagnetics::WireAdviser wireAdviser;
        wireAdviser.enableFoil(true);
        wireAdviser.enableRectangular(false);
        wireAdviser.enableLitz(false);
        wireAdviser.enableRound(false);
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription, section, current, temperature, numberSections, maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(OpenMagnetics::WireType::FOIL == OpenMagnetics::CoilWrapper::resolve_wire(masMagneticWithCoil).get_type());
    }


    TEST(Test_WireAdviser_Low_Frequency_Few_Turns) {
        numberTurns = 2;
        currentRms = 100;
        currentEffectiveFrequency = 13456;
        setup();
        OpenMagnetics::WireAdviser wireAdviser;
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription, section, current, temperature, numberSections, maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(OpenMagnetics::WireType::RECTANGULAR == OpenMagnetics::CoilWrapper::resolve_wire(masMagneticWithCoil).get_type());
    }

    TEST(Test_WireAdviser_Low_Frequency_Many_Turns) {
        numberTurns = 42;
        currentRms = 2;
        currentEffectiveFrequency = 13456;
        setup();
        OpenMagnetics::WireAdviser wireAdviser;
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription,
                                                                 section,
                                                                 current,
                                                                 temperature,
                                                                 numberSections,
                                                                 maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(OpenMagnetics::WireType::ROUND == OpenMagnetics::CoilWrapper::resolve_wire(masMagneticWithCoil).get_type());
    }

    TEST(Test_WireAdviser_Low_Frequency_Gazillion_Turns) {
        numberTurns = 666;
        currentRms = 0.2;
        currentEffectiveFrequency = 13456;
        setup();
        OpenMagnetics::WireAdviser wireAdviser;
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription,
                                                                 section,
                                                                 current,
                                                                 temperature,
                                                                 numberSections,
                                                                 maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(OpenMagnetics::WireType::ROUND == OpenMagnetics::CoilWrapper::resolve_wire(masMagneticWithCoil).get_type());

    }

    TEST(Test_WireAdviser_Medium_Frequency_Few_Turns) {
        numberTurns = 2;
        currentRms = 2;
        currentEffectiveFrequency = 213456;
        setup();
        OpenMagnetics::WireAdviser wireAdviser;
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription,
                                                                 section,
                                                                 current,
                                                                 temperature,
                                                                 numberSections,
                                                                 maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(OpenMagnetics::WireType::LITZ == OpenMagnetics::CoilWrapper::resolve_wire(masMagneticWithCoil).get_type());
        CHECK(OpenMagnetics::CoilWrapper::resolve_wire(masMagneticWithCoil).get_number_conductors().value() < 500);

    }

    TEST(Test_WireAdviser_Medium_High_Frequency_Many_Turns) {
        numberTurns = 42;
        currentRms = 2;
        currentEffectiveFrequency = 613456;
        setup();
        OpenMagnetics::WireAdviser wireAdviser;
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription,
                                                                 section,
                                                                 current,
                                                                 temperature,
                                                                 numberSections,
                                                                 maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(OpenMagnetics::WireType::LITZ == OpenMagnetics::CoilWrapper::resolve_wire(masMagneticWithCoil).get_type());
        CHECK(OpenMagnetics::CoilWrapper::resolve_wire(masMagneticWithCoil).get_number_conductors().value() < 100);
    }

    TEST(Test_WireAdviser_High_Frequency_Few_Turns_High_Current) {
        numberTurns = 5;
        currentRms = 50;
        currentEffectiveFrequency = 4613456;
        setup();
        OpenMagnetics::WireAdviser wireAdviser;
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription,
                                                                 section,
                                                                 current,
                                                                 temperature,
                                                                 numberSections,
                                                                 maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(OpenMagnetics::WireType::LITZ == OpenMagnetics::CoilWrapper::resolve_wire(masMagneticWithCoil).get_type());
    }

    TEST(Test_WireAdviser_High_Frequency_High_Current) {
        numberTurns = 10;
        currentRms = 50;
        currentEffectiveFrequency = 1613456;
        setup();
        OpenMagnetics::WireAdviser wireAdviser;
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription,
                                                                 section,
                                                                 current,
                                                                 temperature,
                                                                 numberSections,
                                                                 maximumNumberResults);
        CHECK(masMagneticsWithCoil.size() > 0);
        if (masMagneticsWithCoil.size() > 0) {
            auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

            CHECK(OpenMagnetics::WireType::RECTANGULAR == OpenMagnetics::CoilWrapper::resolve_wire(masMagneticWithCoil).get_type());
        }
    }
}