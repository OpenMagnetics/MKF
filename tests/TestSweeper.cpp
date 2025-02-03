#include "Sweeper.h"
#include "Painter.h"
#include "Settings.h"
#include "TestingUtils.h"
#include <UnitTest++.h>


SUITE(Sweeper) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");

    TEST(Test_Sweeper_Impedance_Over_Frequency_Many_Turns) {
        settings->set_coil_wind_even_if_not_fit(true);
        std::vector<int64_t> numberTurns = {110, 110};
        std::vector<int64_t> numberParallels = {1, 1};
        std::string shapeName = "T 12.5/7.5/5";
        std::vector<OpenMagnetics::WireWrapper> wires;
        auto wire = OpenMagnetics::find_wire_by_name("Round 0.15 - Grade 1");
        wires.push_back(wire);
        wires.push_back(wire);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
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
        std::vector<OpenMagnetics::CoreGap> gapping;
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = OpenMagnetics::Sweeper().sweep_impedance_over_frequency(magnetic, 1000, 400000, 1000);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Impedance_Over_Frequency_Many_Turns.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, false, true);
        painter.paint_curve(impedanceSweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Impedance_Over_Frequency_Few_Turns) {
        settings->set_coil_wind_even_if_not_fit(true);
        std::vector<int64_t> numberTurns = {18, 18};
        std::vector<int64_t> numberParallels = {1, 1};
        std::string shapeName = "T 12.5/7.5/5";
        std::vector<OpenMagnetics::WireWrapper> wires;
        auto wire = OpenMagnetics::find_wire_by_name("Round 0.425 - Grade 1");
        wires.push_back(wire);
        wires.push_back(wire);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
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
        std::vector<OpenMagnetics::CoreGap> gapping;
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = OpenMagnetics::Sweeper().sweep_impedance_over_frequency(magnetic, 1000, 4000000, 10000);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Impedance_Over_Frequency_Few_Turns.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, false, true);
        painter.paint_curve(impedanceSweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Impedance_Over_Frequency_Larger_Core_Few_Turns) {
        settings->set_coil_wind_even_if_not_fit(true);
        std::vector<int64_t> numberTurns = {9, 9};
        std::vector<int64_t> numberParallels = {1, 1};
        std::string shapeName = "T 36/23/15";
        std::vector<OpenMagnetics::WireWrapper> wires;
        auto wire = OpenMagnetics::find_wire_by_name("Round 2.50 - Grade 1");
        wires.push_back(wire);
        wires.push_back(wire);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
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
        std::vector<OpenMagnetics::CoreGap> gapping;
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = OpenMagnetics::Sweeper().sweep_impedance_over_frequency(magnetic, 1000, 4000000, 100);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Impedance_Over_Frequency_Larger_Core_Few_Turns.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, false, true);
        painter.paint_curve(impedanceSweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Impedance_Over_Frequency_Larger_Core_Many_Turns) {
        settings->set_coil_wind_even_if_not_fit(true);
        std::vector<int64_t> numberTurns = {17, 17};
        std::vector<int64_t> numberParallels = {1, 1};
        std::string shapeName = "T 36/23/15";
        std::vector<OpenMagnetics::WireWrapper> wires;
        auto wire = OpenMagnetics::find_wire_by_name("Round 1.40 - Grade 1");
        wires.push_back(wire);
        wires.push_back(wire);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
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
        std::vector<OpenMagnetics::CoreGap> gapping;
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = OpenMagnetics::Sweeper().sweep_impedance_over_frequency(magnetic, 1000, 4000000, 100);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Impedance_Over_Frequency_Larger_Core_Many_Turns.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, false, true);
        painter.paint_curve(impedanceSweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Resistance_Over_Frequency_Many_Turns) {
        double temperature = 20;
        std::vector<int64_t> numberTurns = {80, 8, 6};
        std::vector<int64_t> numberParallels = {1, 2, 6};
        std::vector<double> turnsRatios = {16, 13};
        std::string shapeName = "ER 28";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        {
            OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("Round T21A01TXXX-1");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         true);

        coil.wind({0, 1, 2}, 1);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C95";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0000008);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = OpenMagnetics::Sweeper().sweep_resistance_over_frequency(magnetic, 0.1, 1000000, 100, 0);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Resistance_Over_Frequency_Many_Turns.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, false, true);
        painter.paint_curve(impedanceSweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Core_Resistance_Over_Frequency_Many_Turns) {
        double temperature = 20;
        std::vector<int64_t> numberTurns = {80, 8, 6};
        std::vector<int64_t> numberParallels = {1, 2, 6};
        std::vector<double> turnsRatios = {16, 13};
        std::string shapeName = "ER 28";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        {
            OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("Round T21A01TXXX-1");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         true);

        coil.wind({0, 1, 2}, 1);

        int64_t numberStacks = 1;
        std::string coreMaterial = "N87";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0000008);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto coreSweep = OpenMagnetics::Sweeper().sweep_core_resistance_over_frequency(magnetic, 10000, 1200000, 20);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Core_Resistance_Over_Frequency_Many_Turns.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, false, true);
        painter.paint_curve(coreSweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Resistance_Over_Frequency_Web_0) {
        OpenMagnetics::MagneticWrapper magnetic = json::parse(R"({"coil":{"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.00635,"columnShape":"rectangular","columnThickness":0,"columnWidth":0.0079375,"coordinates":[0,0,0],"pins":null,"wallThickness":0,"windingWindows":[{"angle":360,"area":0.0007917304360898403,"coordinates":[0.015875,0,0],"height":null,"radialHeight":0.015875,"sectionsAlignment":"spread","sectionsOrientation":"contiguous","shape":"round","width":null}]}},"functionalDescription":[{"connections":null,"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":4,"wire":{"coating":{"breakdownVoltage":1000,"grade":2,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":null,"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.000505,"minimum":0.000495,"nominal":0.0005},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Elektrisola","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 0.5 - Grade 1","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.000544,"minimum":0.0005239999999990001,"nominal":null},"outerHeight":null,"outerWidth":null,"standard":"IEC 60317","standardName":"0.5 mm","strand":null,"type":"round"}},{"connections":null,"isolationSide":"secondary","name":"Secondary","numberParallels":1,"numberTurns":8,"wire":{"coating":{"breakdownVoltage":5000,"grade":2,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":null,"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.001151},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Nearson","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 17.0 - Single Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.001191},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"17 AWG","strand":null,"type":"round"}}],"layersDescription":[{"additionalCoordinates":null,"coordinateSystem":"polar","coordinates":[0.00026700000000000074,89.97744260415998],"dimensions":[0.0005339999999995018,179.95488520831995],"fillingFactor":0.04433090220254894,"insulationMaterial":null,"name":"Primary section 0 layer 0","orientation":"overlapping","partialWindings":[{"connections":null,"parallelsProportion":[1],"winding":"Primary"}],"section":"Primary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveTurns"},{"additionalCoordinates":null,"coordinateSystem":"polar","coordinates":[0.0079375,179.909770417],"dimensions":[0.015875,0.09022958336009401],"fillingFactor":1,"insulationMaterial":null,"name":"Insulation between Primary and Primary section 1 layer 0","orientation":"overlapping","partialWindings":[],"section":"Insulation between Primary and Primary section 1","turnsAlignment":"spread","type":"insulation","windingStyle":null},{"additionalCoordinates":null,"coordinateSystem":"polar","coordinates":[0.0005955000000000005,270.02255739584],"dimensions":[0.0011910000000000035,179.95488520831995],"fillingFactor":0.20659336246431567,"insulationMaterial":null,"name":"Secondary section 0 layer 0","orientation":"overlapping","partialWindings":[{"connections":null,"parallelsProportion":[1],"winding":"Secondary"}],"section":"Secondary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveTurns"},{"additionalCoordinates":null,"coordinateSystem":"polar","coordinates":[0.0079375,359.954885208],"dimensions":[0.015875,0.09022958336009401],"fillingFactor":1,"insulationMaterial":null,"name":"Insulation between Secondary and Secondary section 3 layer 0","orientation":"overlapping","partialWindings":[],"section":"Insulation between Secondary and Secondary section 3","turnsAlignment":"spread","type":"insulation","windingStyle":null}],"sectionsDescription":[{"coordinateSystem":"polar","coordinates":[0.00026700000000000074,89.97744260415998],"dimensions":[0.0005339999999995018,179.95488520831995],"fillingFactor":0.002882066680966784,"layersAlignment":null,"layersOrientation":"overlapping","margin":[0,0],"name":"Primary section 0","partialWindings":[{"connections":null,"parallelsProportion":[1],"winding":"Primary"}],"type":"conduction","windingStyle":"windByConsecutiveTurns"},{"coordinateSystem":"polar","coordinates":[0.0079375,179.9097704166399],"dimensions":[0.015875,0.09022958336009401],"fillingFactor":1,"layersAlignment":null,"layersOrientation":"overlapping","margin":null,"name":"Insulation between Primary and Primary section 1","partialWindings":[],"type":"insulation","windingStyle":null},{"coordinateSystem":"polar","coordinates":[0.0005955000000000005,270.02255739584],"dimensions":[0.0011910000000000035,179.95488520831995],"fillingFactor":0.028673125080251508,"layersAlignment":null,"layersOrientation":"overlapping","margin":[0,0],"name":"Secondary section 0","partialWindings":[{"connections":null,"parallelsProportion":[1],"winding":"Secondary"}],"type":"conduction","windingStyle":"windByConsecutiveTurns"},{"coordinateSystem":"polar","coordinates":[0.0079375,359.95488520832],"dimensions":[0.015875,0.09022958336009401],"fillingFactor":1,"layersAlignment":null,"layersOrientation":"overlapping","margin":null,"name":"Insulation between Secondary and Secondary section 3","partialWindings":[],"type":"insulation","windingStyle":null}],"turnsDescription":[{"additionalCoordinates":[[0.02958105679173332,0.012249463991711979]],"angle":null,"coordinateSystem":"cartesian","coordinates":[0.014420499559776792,0.005971503700616564],"dimensions":[0.0005339999999995,0.0005339999999995],"layer":"Primary section 0 layer 0","length":0.05882761047701695,"name":"Primary parallel 0 turn 0","orientation":"clockwise","parallel":0,"rotation":22.494360651,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":[[0.012261109128712661,0.029576231875169615]],"angle":null,"coordinateSystem":"cartesian","coordinates":[0.00597718060033567,0.014418147456277829],"dimensions":[0.0005339999999995,0.0005339999999995],"layer":"Primary section 0 layer 0","length":0.058827610477016956,"name":"Primary parallel 0 turn 1","orientation":"clockwise","parallel":0,"rotation":67.48308195300001,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":[[-0.012237816956040178,0.029585877123223095]],"angle":null,"coordinateSystem":"cartesian","coordinates":[-0.005965825875312337,0.014422849428093389],"dimensions":[0.0005339999999995,0.0005339999999995],"layer":"Primary section 0 layer 0","length":0.05882761047701695,"name":"Primary parallel 0 turn 2","orientation":"clockwise","parallel":0,"rotation":112.47180325500001,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":[[-0.029571402374279854,0.012272752365237236]],"angle":null,"coordinateSystem":"cartesian","coordinates":[-0.01441579311796108,0.005982856573589742],"dimensions":[0.0005339999999995,0.0005339999999995],"layer":"Primary section 0 layer 0","length":0.05882761047701695,"name":"Primary parallel 0 turn 3","orientation":"clockwise","parallel":0,"rotation":157.460524557,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":[[-0.031719323454203245,-0.006333710591559409]],"angle":null,"coordinateSystem":"cartesian","coordinates":[-0.014983704154163592,-0.0029919441957531033],"dimensions":[0.001191,0.001191],"layer":"Secondary section 0 layer 0","length":0.060891636850425444,"name":"Secondary parallel 0 turn 0","orientation":"clockwise","parallel":0,"rotation":191.292295117,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":[[-0.026882798159381348,-0.01798739929428273]],"angle":null,"coordinateSystem":"cartesian","coordinates":[-0.01269900649166862,-0.00849696147893812],"dimensions":[0.001191,0.001191],"layer":"Secondary section 0 layer 0","length":0.06089163685042544,"name":"Secondary parallel 0 turn 1","orientation":"clockwise","parallel":0,"rotation":213.786655768,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":[[-0.01795563541598219,-0.02690402436547971]],"angle":null,"coordinateSystem":"cartesian","coordinates":[-0.008481956727782842,-0.012709033413994134],"dimensions":[0.001191,0.001191],"layer":"Secondary section 0 layer 0","length":0.060891636850425444,"name":"Secondary parallel 0 turn 2","orientation":"clockwise","parallel":0,"rotation":236.281016419,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":[[-0.006296242411054827,-0.031726782089446674]],"angle":null,"coordinateSystem":"cartesian","coordinates":[-0.0029742448229185584,-0.01498722749488184],"dimensions":[0.001191,0.001191],"layer":"Secondary section 0 layer 0","length":0.060891636850425444,"name":"Secondary parallel 0 turn 3","orientation":"clockwise","parallel":0,"rotation":258.77537707,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":[[0.006321222177019618,-0.031721814583008254]],"angle":null,"coordinateSystem":"cartesian","coordinates":[0.002986044867254216,-0.014984880923809327],"dimensions":[0.001191,0.001191],"layer":"Secondary section 0 layer 0","length":0.06089163685042546,"name":"Secondary parallel 0 turn 4","orientation":"clockwise","parallel":0,"rotation":281.269737721,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":[[0.01797681412055799,-0.02688987772982441]],"angle":null,"coordinateSystem":"cartesian","coordinates":[0.008491961211144235,-0.012702350768201205],"dimensions":[0.001191,0.001191],"layer":"Secondary section 0 layer 0","length":0.060891636850425444,"name":"Secondary parallel 0 turn 5","orientation":"clockwise","parallel":0,"rotation":303.76409837200003,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":[[0.026896953132013648,-0.017966226160890365]],"angle":null,"coordinateSystem":"cartesian","coordinates":[0.012705693075716947,-0.008486959627315216],"dimensions":[0.001191,0.001191],"layer":"Secondary section 0 layer 0","length":0.060891636850425444,"name":"Secondary parallel 0 turn 6","orientation":"clockwise","parallel":0,"rotation":326.25845902300006,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":[[0.031724300794811054,-0.00630873278324186]],"angle":null,"coordinateSystem":"cartesian","coordinates":[0.014986055370741387,-0.0029801450761788814],"dimensions":[0.001191,0.001191],"layer":"Secondary section 0 layer 0","length":0.060891636850425444,"name":"Secondary parallel 0 turn 7","orientation":"clockwise","parallel":0,"rotation":348.7528196740001,"section":"Secondary section 0","winding":"Secondary"}]},"core":{"distributorsInfo":null,"functionalDescription":{"coating":null,"gapping":[],"material":"A07","numberStacks":1,"shape":{"aliases":[],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0635},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.03175},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0127}},"family":"t","familySubtype":null,"magneticCircuit":"closed","name":"T 64/32/12.7","type":"standard"},"type":"toroidal","magneticCircuit":"closed"},"geometricalDescription":[{"coordinates":[0,0,0],"dimensions":null,"insulationMaterial":null,"machining":null,"material":"A07","rotation":[1.5707963267948966,1.5707963267948966,0],"shape":{"aliases":[],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0635},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.03175},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0127}},"family":"t","familySubtype":null,"magneticCircuit":"closed","name":"T 64/32/12.7","type":"standard"},"type":"toroidal"}],"manufacturerInfo":null,"name":"Custom","processedDescription":{"columns":[{"area":0.000202,"coordinates":[0,0,0],"depth":0.0127,"height":0.1496183501272139,"minimumDepth":null,"minimumWidth":null,"shape":"rectangular","type":"central","width":0.015875}],"depth":0.0127,"effectiveParameters":{"effectiveArea":0.00020161250000000003,"effectiveLength":0.14961835012721392,"effectiveVolume":0.000030164929615022918,"minimumArea":0.0002016125},"height":0.0635,"width":0.0635,"windingWindows":[{"angle":360,"area":0.0007917304360898403,"coordinates":[0.015875,0],"height":null,"radialHeight":0.015875,"sectionsAlignment":null,"sectionsOrientation":null,"shape":null,"width":null}]}},"manufacturerInfo":{"name":"OpenMagnetics","reference":"My custom magnetic"}})");

        auto sweep = OpenMagnetics::Sweeper().sweep_resistance_over_frequency(magnetic, 1000, 4e+06, 1000, 0, 25);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Resistance_Over_Frequency_Web_0.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, false, true);
        painter.paint_curve(sweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }
}

