#include "support/Painter.h"
#include "advisers/CoreAdviser.h"
#include "support/Settings.h"
#include "advisers/WireAdviser.h"
#include "processors/Inputs.h"
#include "TestingUtils.h"

#include <UnitTest++.h>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

SUITE(WireAdviser) {
    auto settings = Settings::GetInstance();
    OpenMagnetics::CoilFunctionalDescription coilFunctionalDescription;
    Section section;
    SignalDescriptor current;

    int64_t numberTurns = 42;
    double windingWindowWidth = 0.005;
    double windingWindowHeight = 0.015;
    double currentRms = 2;
    double currentEffectiveFrequency = 13456;
    double temperature = 22;
    uint8_t numberSections = 1;
    size_t maximumNumberResults = 1;


    void setup() {
        coilFunctionalDescription.set_isolation_side(IsolationSide::PRIMARY);
        coilFunctionalDescription.set_name("primary");
        coilFunctionalDescription.set_number_parallels(1);
        coilFunctionalDescription.set_number_turns(numberTurns);
        coilFunctionalDescription.set_wire("Dummy");

        section.set_dimensions({windingWindowWidth, windingWindowHeight});
        section.set_coordinate_system(CoordinateSystem::CARTESIAN);

        Processed processed;
        processed.set_rms(currentRms);
        processed.set_effective_frequency(currentEffectiveFrequency);
        current.set_processed(processed);
    }

    TEST(Test_Round) {
        settings->reset();
        clear_databases();
        numberTurns = 2;
        currentRms = 10;
        currentEffectiveFrequency = 134567;
        setup();
        WireAdviser wireAdviser;
        settings->set_wire_adviser_include_foil(false);
        settings->set_wire_adviser_include_rectangular(false);
        settings->set_wire_adviser_include_litz(false);
        settings->set_wire_adviser_include_round(true);
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription, section, current, temperature, numberSections, maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(masMagneticsWithCoil.size() > 0);

        CHECK(WireType::ROUND == OpenMagnetics::Coil::resolve_wire(masMagneticWithCoil).get_type());
    }


    TEST(Test_Round_IEC_60317) {
        settings->reset();
        clear_databases();
        numberTurns = 2;
        currentRms = 10;
        currentEffectiveFrequency = 134567;
        setup();
        WireAdviser wireAdviser;
        wireAdviser.set_common_wire_standard(WireStandard::IEC_60317);
        settings->set_wire_adviser_include_foil(false);
        settings->set_wire_adviser_include_rectangular(false);
        settings->set_wire_adviser_include_litz(false);
        settings->set_wire_adviser_include_round(true);
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription, section, current, temperature, numberSections, maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(OpenMagnetics::Coil::resolve_wire(masMagneticWithCoil).get_standard().value() == WireStandard::IEC_60317);
        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(WireType::ROUND == OpenMagnetics::Coil::resolve_wire(masMagneticWithCoil).get_type());
    }

    TEST(Test_Round_NEMA_MW_1000_C) {
        settings->reset();
        clear_databases();
        numberTurns = 2;
        currentRms = 10;
        currentEffectiveFrequency = 134567;
        setup();
        WireAdviser wireAdviser;
        wireAdviser.set_common_wire_standard(WireStandard::NEMA_MW_1000_C);
        settings->set_wire_adviser_include_foil(false);
        settings->set_wire_adviser_include_rectangular(false);
        settings->set_wire_adviser_include_litz(false);
        settings->set_wire_adviser_include_round(true);
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription, section, current, temperature, numberSections, maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(OpenMagnetics::Coil::resolve_wire(masMagneticWithCoil).get_standard().value() == WireStandard::NEMA_MW_1000_C);
        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(WireType::ROUND == OpenMagnetics::Coil::resolve_wire(masMagneticWithCoil).get_type());
    }

    TEST(Test_Litz) {
        settings->reset();
        clear_databases();
        numberTurns = 2;
        currentRms = 10;
        currentEffectiveFrequency = 134567;
        setup();
        WireAdviser wireAdviser;
        settings->set_wire_adviser_include_foil(false);
        settings->set_wire_adviser_include_rectangular(false);
        settings->set_wire_adviser_include_litz(true);
        settings->set_wire_adviser_include_round(false);
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription, section, current, temperature, numberSections, maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(WireType::LITZ == OpenMagnetics::Coil::resolve_wire(masMagneticWithCoil).get_type());
    }

    TEST(Test_Rectangular) {
        settings->reset();
        clear_databases();
        numberTurns = 2;
        currentRms = 10;
        currentEffectiveFrequency = 134567;
        setup();
        WireAdviser wireAdviser;
        settings->set_wire_adviser_include_foil(false);
        settings->set_wire_adviser_include_rectangular(true);
        settings->set_wire_adviser_include_litz(false);
        settings->set_wire_adviser_include_round(false);
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription, section, current, temperature, numberSections, maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(WireType::RECTANGULAR == OpenMagnetics::Coil::resolve_wire(masMagneticWithCoil).get_type());
    }

    TEST(Test_Foil) {
        settings->reset();
        clear_databases();
        numberTurns = 2;
        currentRms = 10;
        currentEffectiveFrequency = 134567;
        setup();
        WireAdviser wireAdviser;
        settings->set_wire_adviser_include_foil(true);
        settings->set_wire_adviser_include_rectangular(false);
        settings->set_wire_adviser_include_litz(false);
        settings->set_wire_adviser_include_round(false);
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription, section, current, temperature, numberSections, maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(WireType::FOIL == OpenMagnetics::Coil::resolve_wire(masMagneticWithCoil).get_type());
    }


    TEST(Test_WireAdviser_Low_Frequency_Few_Turns) {
        settings->reset();
        clear_databases();
        numberTurns = 2;
        currentRms = 100;
        currentEffectiveFrequency = 13456;
        setup();
        WireAdviser wireAdviser;
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription, section, current, temperature, numberSections, maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(WireType::RECTANGULAR == OpenMagnetics::Coil::resolve_wire(masMagneticWithCoil).get_type());
    }

    TEST(Test_WireAdviser_Low_Frequency_Many_Turns) {
        settings->reset();
        clear_databases();
        numberTurns = 42;
        currentRms = 2;
        currentEffectiveFrequency = 13456;
        setup();
        WireAdviser wireAdviser;
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription,
                                                                 section,
                                                                 current,
                                                                 temperature,
                                                                 numberSections,
                                                                 maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(WireType::ROUND == OpenMagnetics::Coil::resolve_wire(masMagneticWithCoil).get_type());
    }

    TEST(Test_WireAdviser_Low_Frequency_Gazillion_Turns) {
        settings->reset();
        clear_databases();
        numberTurns = 666;
        currentRms = 0.2;
        currentEffectiveFrequency = 13456;
        setup();
        WireAdviser wireAdviser;
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription,
                                                                 section,
                                                                 current,
                                                                 temperature,
                                                                 numberSections,
                                                                 maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(WireType::ROUND == OpenMagnetics::Coil::resolve_wire(masMagneticWithCoil).get_type());

    }

    TEST(Test_WireAdviser_Medium_Frequency_Few_Turns) {
        settings->reset();
        clear_databases();
        numberTurns = 2;
        currentRms = 2;
        currentEffectiveFrequency = 213456;
        setup();
        WireAdviser wireAdviser;
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription,
                                                                 section,
                                                                 current,
                                                                 temperature,
                                                                 numberSections,
                                                                 maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(WireType::LITZ == OpenMagnetics::Coil::resolve_wire(masMagneticWithCoil).get_type());
        CHECK(OpenMagnetics::Coil::resolve_wire(masMagneticWithCoil).get_number_conductors().value() < 500);

    }

    TEST(Test_WireAdviser_Medium_High_Frequency_Many_Turns) {
        settings->reset();
        clear_databases();
        numberTurns = 42;
        currentRms = 2;
        currentEffectiveFrequency = 613456;
        setup();
        WireAdviser wireAdviser;
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription,
                                                                 section,
                                                                 current,
                                                                 temperature,
                                                                 numberSections,
                                                                 maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(WireType::LITZ == OpenMagnetics::Coil::resolve_wire(masMagneticWithCoil).get_type());
        CHECK(OpenMagnetics::Coil::resolve_wire(masMagneticWithCoil).get_number_conductors().value() < 100);
    }

    TEST(Test_WireAdviser_High_Frequency_Few_Turns_High_Current) {
        settings->reset();
        clear_databases();
        numberTurns = 5;
        currentRms = 50;
        currentEffectiveFrequency = 4613456;
        setup();
        WireAdviser wireAdviser;
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription,
                                                                 section,
                                                                 current,
                                                                 temperature,
                                                                 numberSections,
                                                                 maximumNumberResults);
        auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

        CHECK(masMagneticsWithCoil.size() > 0);
        CHECK(WireType::LITZ == OpenMagnetics::Coil::resolve_wire(masMagneticWithCoil).get_type());
    }

    TEST(Test_WireAdviser_High_Frequency_High_Current) {
        settings->reset();
        clear_databases();
        numberTurns = 10;
        currentRms = 50;
        currentEffectiveFrequency = 1613456;
        setup();
        WireAdviser wireAdviser;
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription,
                                                                 section,
                                                                 current,
                                                                 temperature,
                                                                 numberSections,
                                                                 maximumNumberResults);
        CHECK(masMagneticsWithCoil.size() > 0);
        if (masMagneticsWithCoil.size() > 0) {
            auto masMagneticWithCoil = masMagneticsWithCoil[0].first;

            CHECK(WireType::RECTANGULAR == OpenMagnetics::Coil::resolve_wire(masMagneticWithCoil).get_type());
        }
    }

    TEST(Test_WireAdviser_Web_0) {
        settings->reset();
        clear_databases();
        settings->set_painter_simple_litz(false);

        WireAdviser wireAdviser;
        OpenMagnetics::CoilFunctionalDescription coilFunctionalDescription(json::parse(R"({"connections":null,"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":1,"wire":"Dummy"})"));
        Section section(json::parse(R"({"coordinateSystem":"cartesian","coordinates":[0.010311213441920776,0,0],"dimensions":[0.001395915,0.017288338192419827],"fillingFactor":0.0000064745233637820475,"layersAlignment":null,"layersOrientation":"overlapping","margin":[0,0],"name":"Primary section 0","partialWindings":[{"connections":null,"parallelsProportion":[1],"winding":"Primary"}],"type":"conduction","windingStyle":"windByConsecutiveParallels"})"));
        SignalDescriptor current(json::parse(R"({"waveform":{"ancillaryLabel":null,"data":[-5,5,-5],"numberPeriods":null,"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular","acEffectiveFrequency":110746.40291779551,"effectiveFrequency":110746.40291779551,"peak":5,"rms":2.8874560332150576,"thd":0.12151487440704967},"harmonics":{"amplitudes":[1.1608769501236793e-14,4.05366124583194,1.787369544444173e-15,0.4511310569983995,9.749053004706756e-16,0.16293015292554872,4.036157626725542e-16,0.08352979924600704,3.4998295008010614e-16,0.0508569581336163,3.1489164048780735e-16,0.034320410449418075,3.142469873118059e-16,0.024811988673843106,2.3653352035940994e-16,0.018849001010678823,2.9306524147249266e-16,0.014866633059596499,1.796485796132569e-16,0.012077180559557785,1.6247782523152451e-16,0.010049063750920609,1.5324769149805092e-16,0.008529750975091871,1.0558579733068502e-16,0.007363501410705499,7.513269775674661e-17,0.006450045785294609,5.871414177162291e-17,0.005722473794997712,9.294731722001391e-17,0.005134763398167541,1.194820309200107e-16,0.004654430423785411,8.2422739080512e-17,0.004258029771397032,9.5067306351894e-17,0.0039283108282380024,1.7540347128474968e-16,0.0036523670873925395,9.623794010508822e-17,0.0034204021424253787,1.4083470894369491e-16,0.0032248884817922927,1.4749333016985644e-16,0.0030599828465501895,1.0448590642474364e-16,0.002921112944200383,7.575487373767413e-18,0.002804680975178716,7.419510610361002e-17,0.0027078483284668376,3.924741709073613e-17,0.0026283777262804953,2.2684279102637236e-17,0.0025645167846443107,8.997077625295079e-17,0.0025149120164513483,7.131074184849219e-17,0.0024785457043284276,9.354417496250849e-17,0.0024546904085875065,1.2488589642405877e-16,0.0024428775264784264],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}})"));
        auto masMagneticsWithCoil = wireAdviser.get_advised_wire(coilFunctionalDescription,
                                                                 section,
                                                                 current,
                                                                 25,
                                                                 1,
                                                                 1);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_WireAdviser_Web_0.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            auto wire = OpenMagnetics::Coil::resolve_wire(masMagneticsWithCoil[0].first);
            InsulationWireCoating newCoating;
            newCoating.set_type(InsulationWireCoatingType::INSULATED);
            newCoating.set_number_layers(2);
            newCoating.set_thickness_layers(5.08e-05);
            wire.set_coating(newCoating);

            painter.paint_wire(wire);
            painter.export_svg();
        }


    }
}