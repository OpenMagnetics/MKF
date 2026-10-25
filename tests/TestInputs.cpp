#include <source_location>
#include "support/CoilMesher.h"
#include "support/Painter.h"
#include "physical_models/InitialPermeability.h"
#include "processors/Inputs.h"
#include "processors/WaveformProcessor.h"
#include "processors/ResultsAnalyzer.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/NumberFormat.h"
#include <regex>
#include "support/Utils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
using json = nlohmann::json;
#include "TestingUtils.h"
#include "Fixtures.h"

#include <typeinfo>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace MAS;
using namespace OpenMagnetics;

namespace {
auto outputFilePath = std::filesystem::path {std::source_location::current().file_name()}.parent_path().append("..").append("output");

TEST_CASE("Test_One_Operating_Point_One_Winding_Triangular", "[processor][inputs][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["current"]["waveform"]["data"] = {-5, 5, -5};
    windingExcitation["current"]["waveform"]["time"] = {0, 0.0000025, 0.00001};
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);
    double max_error = 0.01;

    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_rms().value(), Catch::Matchers::WithinAbs(2.88, max_error * 2.88));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_thd().value(), Catch::Matchers::WithinAbs(0.376, max_error * 0.376));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(128000, max_error * 128000));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_peak_to_peak().value(), Catch::Matchers::WithinAbs(10, max_error * 10));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_offset(), Catch::Matchers::WithinAbs(0, max_error));

    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_amplitudes()[0], Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_amplitudes()[1], Catch::Matchers::WithinAbs(3.833, max_error * 3.833));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_frequencies()[0], Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_frequencies()[1], Catch::Matchers::WithinAbs(100000, max_error * 100000));
}

TEST_CASE("Test_One_Operating_Point_One_Winding_Equidistant_Sinusoidal", "[processor][inputs][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["current"]["waveform"]["data"] = {0, 0.06145122063750744, 0.12290244127501487, 0.18435366191252228, 0.24580488255002975, 0.30713118165698094, 0.36841396877014293, 0.4296967558833049, 0.49097954299646696, 0.5520997341739765, 0.6130461159030302, 0.6739924976320839, 0.7349388793611376, 0.7957719836746695, 0.8562149102236266, 0.9166578367725836, 0.9771007633215409, 1.037543689870498, 1.0973397955470385, 1.1571135970555646, 1.216887398564091, 1.276661200072617, 1.3358449876276357, 1.394785828260622, 1.4537266688936081, 1.5126675095265945, 1.5711608190293265, 1.6291071460395254, 1.6870534730497244, 1.7449998000599234, 1.8027250701453021, 1.8595180566781788, 1.9163110432110557, 1.973104029743932, 2.0298970162768093, 2.085468263777627, 2.140952244204727, 2.1964362246318263, 2.251920205058926, 2.3064179538943743, 2.3604408504819907, 2.4144637470696066, 2.468486643657223, 2.5218390582265293, 2.574252797970201, 2.6266665377138723, 2.6790802774575444, 2.731216487484836, 2.78187740796716, 2.8325383284494845, 2.883199248931809, 2.933860169414133, 2.9828185803057643, 3.031587823465096, 3.0803570666244284, 3.1291263097837603, 3.176595953084751, 3.223339845815499, 3.2700837385462473, 3.3168276312769955, 3.362746026335129, 3.4073364468764473, 3.4519268674177654, 3.496517287959084, 3.540823247071316, 3.5831379761799416, 3.625452705288568, 3.6677674343971947, 3.710082163505821, 3.750324109205718, 3.790247165150427, 3.8301702210951354, 3.8700932770398437, 3.9084948307976948, 3.9459167872571586, 3.983338743716622, 4.020760700176085, 4.057271371972193, 4.092089657964925, 4.126907943957656, 4.161726229950387, 4.196297097987679, 4.228416279011962, 4.2605354600362455, 4.292654641060529, 4.324773822084812, 4.3545703852665465, 4.38390242487977, 4.413234464492994, 4.442566504106219, 4.470249709086087, 4.496714210206516, 4.523178711326947, 4.549643212447377, 4.575176689548374, 4.598701114817964, 4.6222255400875545, 4.645749965357144, 4.669099124937923, 4.689618995542171, 4.71013886614642, 4.730658736750669, 4.751178607354918, 4.7692498393956875, 4.7867089117955866, 4.8041679841954865, 4.8216270565953865, 4.8374022756616615, 4.85175269575333, 4.866103115844998, 4.880453535936666, 4.893912026714872, 4.905114460994323, 4.9163168952737735, 4.927519329553225, 4.938642296560028, 4.9466660399335725, 4.954689783307117, 4.962713526680662, 4.970737270054206, 4.976307156132361, 4.981130216077137, 4.985953276021912, 4.990776335966688, 4.993965661835207, 4.9955748186791515, 4.997183975523097, 4.998793132367041, 4.999597710789014, 4.99798855394507, 4.9963793971011246, 4.99477024025718, 4.993161083413235, 4.9883648059942995, 4.983541746049524, 4.978718686104749, 4.973895626159973, 4.9667253983674335, 4.95870165499389, 4.950677911620344, 4.942654168246801, 4.93312054669295, 4.9219181124135, 4.910715678134049, 4.899513243854598, 4.8876287459825, 4.873278325890832, 4.858927905799163, 4.844577485707495, 4.830227065615827, 4.812897520395438, 4.795438447995538, 4.777979375595638, 4.760520303195739, 4.7409186720527945, 4.720398801448546, 4.699878930844296, 4.679359060240047, 4.657512177991942, 4.633987752722351, 4.610463327452761, 4.586938902183172, 4.562875463007594, 4.536410961887164, 4.5099464607667334, 4.483481959646303, 4.457017458525873, 4.427900484299609, 4.398568444686385, 4.36923640507316, 4.339904365459936, 4.308714231572671, 4.276595050548388, 4.244475869524105, 4.212356688499822, 4.1791353729467575, 4.144317086954025, 4.109498800961293, 4.074680514968562, 4.039471678405818, 4.002049721946355, 3.96462776548689, 3.927205809027427, 3.8897838525679624, 3.85013174906749, 3.8102086931227817, 3.7702856371780733, 3.730362581233366, 3.6889247989515104, 3.6466100698428843, 3.6042953407342573, 3.5619806116256303, 3.518812498229744, 3.4742220776884256, 3.429631657147107, 3.3850412366057903, 3.3401995776423714, 3.2934556849116223, 3.246711792180875, 3.199967899450126, 3.1532240067193786, 3.1047416882040952, 3.0559724450447625, 3.0072032018854333, 2.9584339587261006, 2.9085297091729725, 2.8578687886906486, 2.807207868208325, 2.756546947726001, 2.7052871473293862, 2.652873407585714, 2.60045966784204, 2.548045928098368, 2.495498091951033, 2.44147519536342, 2.387452298775802, 2.3334294021881856, 2.2794065056005692, 2.2241782148453773, 2.1686942344182754, 2.113210253991177, 2.057726273564077, 2.001500523010371, 1.9447075364774946, 1.88791454994462, 1.8311215634117453, 1.7739729635650274, 1.7160266365548278, 1.6580803095446264, 1.6001339825344267, 1.5421379298430882, 1.483197089210103, 1.4242562485771177, 1.3653154079441308, 1.3063745673111438, 1.2467742993183535, 1.18700049780983, 1.127226696301305, 1.067452894792778, 1.0073222265960222, 0.9468793000470654, 0.8864363734981069, 0.8259934469491519, 0.7654120702256666, 0.7044656884966134, 0.6435193067675584, 0.5825729250385052, 0.5216209365530489, 0.46033814943988993, 0.39905536232672745, 0.33777257521356496, 0.27648978810040425, 0.21507927223127865, 0.15362805159377224, 0.09217683095626406, 0.03072561031875587, -0.03072561031875054, -0.09217683095625517, -0.15362805159376514, -0.21507927223126977, -0.2764897881003954, -0.33777257521355786, -0.39905536232672034, -0.4603381494398846, -0.5216209365530453, -0.5825729250384999, -0.6435193067675531, -0.7044656884966063, -0.765412070225663, -0.8259934469491466, -0.8864363734981016, -0.9468793000470601, -1.0073222265960187, -1.0674528947927726, -1.1272266963012996, -1.187000497809823, -1.24677429931835, -1.3063745673111367, -1.3653154079441254, -1.4242562485771106, -1.4831970892100994, -1.542137929843081, -1.6001339825344232, -1.6580803095446228, -1.7160266365548225, -1.773972963565022, -1.8311215634117364, -1.8879145499446164, -1.9447075364774928, -2.0015005230103693, -2.057726273564077, -2.1132102539911752, -2.1686942344182736, -2.224178214845372, -2.279406505600564, -2.333429402188182, -2.3874522987757985, -2.441475195363413, -2.4954980919510312, -2.5480459280983627, -2.6004596678420366, -2.6528734075857088, -2.705287147329381, -2.7565469477259974, -2.8072078682083212, -2.857868788690647, -2.908529709172969, -2.9584339587260953, -3.007203201885428, -3.055972445044759, -3.1047416882040917, -3.153224006719375, -3.199967899450124, -3.2467117921808732, -3.2934556849116206, -3.3401995776423696, -3.3850412366057885, -3.4296316571471053, -3.474222077688424, -3.5188124982297424, -3.5619806116256267, -3.6042953407342537, -3.6466100698428807, -3.688924798951504, -3.7303625812333614, -3.770285637178068, -3.810208693122778, -3.8501317490674865, -3.8897838525679607, -3.9272058090274253, -3.96462776548689, -4.002049721946351, -4.039471678405816, -4.074680514968556, -4.109498800961289, -4.144317086954022, -4.179135372946753, -4.2123566884998205, -4.244475869524104, -4.276595050548385, -4.3087142315726705, -4.339904365459933, -4.369236405073156, -4.398568444686381, -4.427900484299604, -4.45701745852587, -4.4834819596463, -4.509946460766729, -4.53641096188716, -4.56287546300759, -4.586938902183168, -4.610463327452758, -4.633987752722349, -4.657512177991937, -4.6793590602400466, -4.699878930844295, -4.7203988014485425, -4.740918672052792, -4.760520303195736, -4.777979375595635, -4.795438447995536, -4.812897520395436, -4.830227065615827, -4.844577485707496, -4.858927905799163, -4.873278325890832, -4.887628745982499, -4.899513243854597, -4.910715678134047, -4.921918112413499, -4.933120546692948, -4.942654168246799, -4.950677911620344, -4.958701654993888, -4.9667253983674335, -4.973895626159973, -4.978718686104749, -4.983541746049525, -4.9883648059943, -4.993161083413236, -4.99477024025718, -4.9963793971011246, -4.99798855394507, -4.999597710789014, -4.998793132367042, -4.997183975523097, -4.995574818679152, -4.993965661835208, -4.990776335966688, -4.985953276021913, -4.981130216077137, -4.9763071561323615, -4.970737270054206, -4.962713526680662, -4.954689783307117, -4.9466660399335725, -4.938642296560029, -4.927519329553226, -4.916316895273774, -4.905114460994324, -4.893912026714872, -4.880453535936667, -4.866103115844997, -4.851752695753331, -4.837402275661661, -4.821627056595389, -4.804167984195489, -4.786708911795589, -4.769249839395689, -4.75117860735492, -4.730658736750671, -4.7101388661464245, -4.689618995542174, -4.669099124937926, -4.6457499653571475, -4.622225540087557, -4.598701114817967, -4.5751766895483765, -4.54964321244738, -4.523178711326947, -4.49671421020652, -4.470249709086087, -4.442566504106223, -4.413234464493, -4.383902424879775, -4.35457038526655, -4.324773822084817, -4.2926546410605315, -4.260535460036252, -4.228416279011965, -4.196297097987683, -4.161726229950389, -4.126907943957658, -4.09208965796493, -4.057271371972195, -4.020760700176087, -3.9833387437166223, -3.9459167872571577, -3.908494830797693, -3.870093277039846, -3.8301702210951376, -3.7902471651504257, -3.7503241092057173, -3.7100821635058203, -3.667767434397195, -3.6254527052885734, -3.583137976179941, -3.5408232470713195, -3.4965172879590867, -3.451926867417768, -3.407336446876453, -3.362746026335131, -3.316827631277004, -3.2700837385462513, -3.223339845815506, -3.176595953084753, -3.1291263097837643, -3.0803570666244298, -3.031587823465099, -2.9828185803057643, -2.9338601694141317, -2.8831992489318097, -2.8325383284494876, -2.7818774079671655, -2.73121648748484, -2.6790802774575475, -2.626666537713877, -2.574252797970207, -2.5218390582265293, -2.4684866436572293, -2.4144637470696075, -2.3604408504819965, -2.3064179538943783, -2.25192020505893, -2.196436224631828, -2.1409522442047297, -2.0854682637776314, -2.029897016276813, -1.9731040297439364, -1.9163110432110635, -1.859518056678187, -1.802725070145307, -1.7449998000599258, -1.6870534730497262, -1.62910714603953, -1.571160819029327, -1.512667509526601, -1.4537266688936121, -1.3947858282606305, -1.3358449876276417, -1.276661200072624, -1.2168873985640971, -1.1571135970555666, -1.0973397955470396, -1.0375436898704997, -0.9771007633215412, -0.9166578367725862, -0.8562149102236276, -0.7957719836746691, -0.7349388793611382, -0.6739924976320886, -0.6130461159030389, -0.5520997341739822, -0.4909795429964774, -0.4296967558833096, -0.3684139687701524, -0.30713118165698816, -0.24580488255003274, -0.18435366191252456, -0.12290244127501992, -0.06145122063751174, 0};
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);
    double max_error = 0.01;

    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_rms().value(), Catch::Matchers::WithinAbs(3.53, max_error * 3.53));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_thd().value(), Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(100000, max_error * 100000));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_peak_to_peak().value(), Catch::Matchers::WithinAbs(10, max_error * 10));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_offset(), Catch::Matchers::WithinAbs(0, max_error));

    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_amplitudes()[0], Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_amplitudes()[1], Catch::Matchers::WithinAbs(5, max_error * 5));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_frequencies()[0], Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_frequencies()[1], Catch::Matchers::WithinAbs(100000, max_error * 100000));
}

TEST_CASE("Test_One_Operating_Point_One_Winding_Equidistant_Rectangular", "[processor][inputs][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["voltage"]["waveform"]["data"] = {-2.5, 7.5, 7.5, -2.5, -2.5};
    windingExcitation["voltage"]["waveform"]["time"] = {0, 0, 0.0000025, 0.0000025, 0.00001};
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);
    auto max_error = 0.1;

    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    REQUIRE_THAT(excitation.get_voltage().value().get_processed().value().get_rms().value(), Catch::Matchers::WithinAbs(4.28, max_error * 4.28));
    REQUIRE_THAT(excitation.get_voltage().value().get_processed().value().get_thd().value(), Catch::Matchers::WithinAbs(0.95, max_error * 0.95));
    // ABT #1460: exact 64-harmonic series of the piecewise-linear waveform (previously an aliased DFT128).
    REQUIRE_THAT(excitation.get_voltage().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(590582.20189317234, max_error * 590582.20189317234));
    REQUIRE_THAT(excitation.get_voltage().value().get_processed().value().get_peak_to_peak().value(), Catch::Matchers::WithinAbs(10, max_error * 10));
    REQUIRE_THAT(excitation.get_voltage().value().get_processed().value().get_offset(), Catch::Matchers::WithinAbs(0, max_error));

    REQUIRE_THAT(excitation.get_voltage().value().get_harmonics().value().get_amplitudes()[0], Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_voltage().value().get_harmonics().value().get_amplitudes()[1], Catch::Matchers::WithinAbs(4.5, max_error * 4.5));
    REQUIRE_THAT(excitation.get_voltage().value().get_harmonics().value().get_frequencies()[0], Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_voltage().value().get_harmonics().value().get_frequencies()[1], Catch::Matchers::WithinAbs(100000, max_error * 100000));
}

TEST_CASE("Test_One_Operating_Point_One_Winding_Triangular_Processed", "[processor][inputs][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["current"]["processed"]["dutyCycle"] = 0.8;
    windingExcitation["current"]["processed"]["label"] = WaveformLabel::TRIANGULAR;
    windingExcitation["current"]["processed"]["offset"] = 5;
    windingExcitation["current"]["processed"]["peakToPeak"] = 20;
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);
    double max_error = 0.01;

    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_rms().value(), Catch::Matchers::WithinAbs(7.64, max_error * 7.64));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_thd().value(), Catch::Matchers::WithinAbs(0.452, max_error * 0.452));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(117400, max_error * 117400));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_peak_to_peak().value(), Catch::Matchers::WithinAbs(20, max_error * 20));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_offset(), Catch::Matchers::WithinAbs(5, max_error * 5));

    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_amplitudes()[0], Catch::Matchers::WithinAbs(5, max_error));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_amplitudes()[1], Catch::Matchers::WithinAbs(7.45, max_error * 7.45));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_frequencies()[0], Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_frequencies()[1], Catch::Matchers::WithinAbs(100000, max_error * 100000));
}

TEST_CASE("Test_One_Operating_Point_One_Winding_Rectangular_Processed", "[processor][inputs][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["current"]["processed"]["dutyCycle"] = 0.8;
    windingExcitation["current"]["processed"]["label"] = WaveformLabel::RECTANGULAR;
    windingExcitation["current"]["processed"]["offset"] = 0;
    windingExcitation["current"]["processed"]["peakToPeak"] = 20;
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);
    double max_error = 0.1;

    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_rms().value(), Catch::Matchers::WithinAbs(7.93, max_error * 7.93));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_thd().value(), Catch::Matchers::WithinAbs(1.15, max_error * 1.15));
    // ABT #1460: exact 64-harmonic series of the piecewise-linear waveform (previously an aliased DFT128).
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(644721.39218268555, max_error * 644721.39218268555));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_peak_to_peak().value(), Catch::Matchers::WithinAbs(20, max_error * 20));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_offset(), Catch::Matchers::WithinAbs(0, max_error));

    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_amplitudes()[0], Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_amplitudes()[1], Catch::Matchers::WithinAbs(7.33, max_error * 7.33));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_frequencies()[0], Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_frequencies()[1], Catch::Matchers::WithinAbs(100000, max_error * 100000));
}

TEST_CASE("Test_One_Operating_Point_One_Winding_Sinusoidal_Processed", "[processor][inputs][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["current"]["processed"]["label"] = WaveformLabel::SINUSOIDAL;
    windingExcitation["current"]["processed"]["offset"] = -6;
    windingExcitation["current"]["processed"]["peakToPeak"] = 12;
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);
    double max_error = 0.1;

    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_rms().value(), Catch::Matchers::WithinAbs(7.35, max_error * 7.35));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_thd().value(), Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(70700, max_error * 70700));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_peak_to_peak().value(), Catch::Matchers::WithinAbs(12, max_error * 12));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_offset() + 6, Catch::Matchers::WithinAbs(0, max_error));

    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_amplitudes()[0], Catch::Matchers::WithinAbs(6, max_error * 6));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_amplitudes()[1], Catch::Matchers::WithinAbs(6, max_error * 6));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_frequencies()[0], Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_frequencies()[1], Catch::Matchers::WithinAbs(100000, max_error * 100000));
}

TEST_CASE("Test_One_Operating_Point_One_Winding_Bipolar_Rectangular_Processed", "[processor][inputs][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["current"]["processed"]["dutyCycle"] = 0.42;
    windingExcitation["current"]["processed"]["label"] = WaveformLabel::BIPOLAR_RECTANGULAR;
    windingExcitation["current"]["processed"]["offset"] = 0;
    windingExcitation["current"]["processed"]["peakToPeak"] = 13;
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);
    double max_error = 0.01;

    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_rms().value(), Catch::Matchers::WithinAbs(5.92, max_error * 5.92));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_thd().value(), Catch::Matchers::WithinAbs(0.315, max_error * 0.315));
    // ABT #1460: exact 64-harmonic series of the piecewise-linear waveform (previously an aliased DFT128).
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(398003.45508438692, max_error * 398003.45508438692));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_peak_to_peak().value(), Catch::Matchers::WithinAbs(13, max_error * 13));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_offset(), Catch::Matchers::WithinAbs(0, max_error));

    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_amplitudes()[0], Catch::Matchers::WithinAbs(0.0, max_error * 0.0));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_amplitudes()[1], Catch::Matchers::WithinAbs(8, max_error * 8));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_frequencies()[0], Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_frequencies()[1], Catch::Matchers::WithinAbs(100000, max_error * 100000));
}

TEST_CASE("Test_One_Operating_Point_One_Winding_Bipolar_Triangular_Processed", "[processor][inputs][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["current"]["processed"]["dutyCycle"] = 0.42;
    windingExcitation["current"]["processed"]["label"] = WaveformLabel::BIPOLAR_TRIANGULAR;
    windingExcitation["current"]["processed"]["offset"] = 0;
    windingExcitation["current"]["processed"]["peakToPeak"] = 13;
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);
    double max_error = 0.01;

    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_rms().value(), Catch::Matchers::WithinAbs(4.31, max_error * 4.31));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_thd().value(), Catch::Matchers::WithinAbs(0.086, max_error * 0.086));
    REQUIRE_THAT(104000, Catch::Matchers::WithinAbs(excitation.get_current().value().get_processed().value().get_effective_frequency().value(), max_error * 104000));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_peak_to_peak().value(), Catch::Matchers::WithinAbs(13, max_error * 13));
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_offset(), Catch::Matchers::WithinAbs(0, max_error));

    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_amplitudes()[0], Catch::Matchers::WithinAbs(0.0, max_error));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_amplitudes()[1], Catch::Matchers::WithinAbs(6.07, max_error * 6.07));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_frequencies()[0], Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_current().value().get_harmonics().value().get_frequencies()[1], Catch::Matchers::WithinAbs(100000, max_error * 100000));
}

TEST_CASE("Test_One_Operating_Point_Two_Generated_Windings_Turns_Ratios", "[processor][inputs][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;
    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["current"]["waveform"]["data"] = {-5, 5, -5};
    windingExcitation["current"]["waveform"]["time"] = {0, 0.0000025, 0.00001};
    // Volt-second balanced (ABT #1585): the same 1 V peak-to-peak, D = 0.25 square, shifted to zero mean
    // (+0.75 V for 25 %, -0.25 V for 75 %); the former unipolar 0/1 V square is rejected as unbalanced.
    windingExcitation["voltage"]["waveform"]["data"] = {-0.25, 0.75, 0.75, -0.25, -0.25};
    windingExcitation["voltage"]["waveform"]["time"] = {0, 0, 0.0000025, 0.0000025, 0.00001};
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    std::vector<double> turnsRatioValues = {0.001, 0.1, 1, 5, 42};
    auto max_error = 0.02;

    for (auto& turnsRatioValue : turnsRatioValues) {
        inputsJson["designRequirements"] = json();
        inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
        json turnsRatio;
        turnsRatio["nominal"] = turnsRatioValue;
        inputsJson["designRequirements"]["turnsRatios"] = json::array({turnsRatio});

        OpenMagnetics::Inputs inputs(inputsJson);

        auto excitation_primary = inputs.get_operating_points()[0].get_excitations_per_winding()[0];

        REQUIRE_THAT(excitation_primary.get_current().value().get_processed().value().get_rms().value(), Catch::Matchers::WithinAbs(2.9, max_error * 2.9));
        REQUIRE_THAT(excitation_primary.get_current().value().get_processed().value().get_thd().value(), Catch::Matchers::WithinAbs(0.382, max_error * 0.382));
        // ABT #1460: exact 64-harmonic series of the piecewise-linear waveform (previously an aliased DFT128).
        REQUIRE_THAT(excitation_primary.get_current().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(126785.46967026753, max_error * 126785.46967026753));
        REQUIRE_THAT(excitation_primary.get_current().value().get_processed().value().get_peak_to_peak().value(), Catch::Matchers::WithinAbs(10, max_error * 10));
        REQUIRE_THAT(excitation_primary.get_current().value().get_processed().value().get_offset(), Catch::Matchers::WithinAbs(0, max_error));

        REQUIRE_THAT(excitation_primary.get_current().value().get_harmonics().value().get_amplitudes()[0], Catch::Matchers::WithinAbs(0, max_error));
        REQUIRE_THAT(excitation_primary.get_current().value().get_harmonics().value().get_amplitudes()[1], Catch::Matchers::WithinAbs(3.833, max_error * 3.833));
        REQUIRE_THAT(excitation_primary.get_current().value().get_harmonics().value().get_frequencies()[0], Catch::Matchers::WithinAbs(0, max_error));
        REQUIRE_THAT(excitation_primary.get_current().value().get_harmonics().value().get_frequencies()[1], Catch::Matchers::WithinAbs(100000, max_error * 100000));

        // ABT #1585: the balanced fixture (-0.25/+0.75 V, D = 0.25) has RMS sqrt(0.25 * 0.75^2 + 0.75 * 0.25^2) = 0.433 V
        // (the former unipolar 0/1 V square: 0.5 V). The processed value, 0.428 V, is the 128-point sampled waveform's.
        REQUIRE_THAT(excitation_primary.get_voltage().value().get_processed().value().get_rms().value(), Catch::Matchers::WithinAbs(0.43301270189221935, max_error * 0.433));
        // ABT #1460: exact 64-harmonic series of the piecewise-linear waveform (previously an aliased DFT128).
        REQUIRE_THAT(excitation_primary.get_voltage().value().get_processed().value().get_thd().value(), Catch::Matchers::WithinAbs(0.91374449406562508, max_error * 0.92));
        // ABT #1460: exact 64-harmonic series of the piecewise-linear waveform (previously an aliased DFT128).
        // ABT #1585: the balanced fixture has no DC component, which the unipolar 0/1 V square's 546440.68 Hz included.
        REQUIRE_THAT(excitation_primary.get_voltage().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(590582.20189317246, max_error * 590582.20189317246));
        REQUIRE_THAT(excitation_primary.get_voltage().value().get_processed().value().get_peak_to_peak().value(), Catch::Matchers::WithinAbs(1, max_error * 1));
        REQUIRE_THAT(excitation_primary.get_voltage().value().get_processed().value().get_offset(), Catch::Matchers::WithinAbs(0, max_error * 0.25));

        // ABT #1585: the balanced fixture's mean is 0.25 * 0.75 - 0.75 * 0.25 = 0 V (the former unipolar 0/1 V square's DC was 0.25 V).
        REQUIRE_THAT(excitation_primary.get_voltage().value().get_harmonics().value().get_amplitudes()[0], Catch::Matchers::WithinAbs(0, max_error * 0.25));
        // ABT #1460: exact 64-harmonic series of the piecewise-linear waveform (previously an aliased DFT128). Fundamental = (2/pi) sin(pi/4).
        REQUIRE_THAT(excitation_primary.get_voltage().value().get_harmonics().value().get_amplitudes()[1], Catch::Matchers::WithinAbs(0.45015815807855303, max_error * 0.451));
        REQUIRE_THAT(excitation_primary.get_voltage().value().get_harmonics().value().get_frequencies()[0], Catch::Matchers::WithinAbs(0, max_error));
        REQUIRE_THAT(excitation_primary.get_voltage().value().get_harmonics().value().get_frequencies()[1], Catch::Matchers::WithinAbs(100000, max_error * 100000));

        auto excitation_secondary = inputs.get_operating_points()[0].get_excitations_per_winding()[1];
        REQUIRE_THAT(excitation_secondary.get_current().value().get_processed().value().get_rms().value(), Catch::Matchers::WithinAbs(excitation_primary.get_current().value().get_processed().value().get_rms().value() *
                        turnsRatioValue, max_error * excitation_primary.get_current().value().get_processed().value().get_rms().value() *
                        turnsRatioValue));
        REQUIRE_THAT(excitation_secondary.get_current().value().get_processed().value().get_thd().value(), Catch::Matchers::WithinAbs(excitation_primary.get_current().value().get_processed().value().get_thd().value(), max_error * excitation_primary.get_current().value().get_processed().value().get_thd().value()));
        REQUIRE_THAT(excitation_secondary.get_current().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(excitation_primary.get_current().value().get_processed().value().get_effective_frequency().value(), max_error *
                excitation_primary.get_current().value().get_processed().value().get_effective_frequency().value()));
        REQUIRE_THAT(excitation_secondary.get_current().value().get_processed().value().get_peak_to_peak().value(), Catch::Matchers::WithinAbs(excitation_primary.get_current().value().get_processed().value().get_peak_to_peak().value() *
                turnsRatioValue, max_error *
                excitation_primary.get_current().value().get_processed().value().get_peak_to_peak().value() *
                turnsRatioValue));
        REQUIRE_THAT(excitation_secondary.get_current().value().get_processed().value().get_offset(), Catch::Matchers::WithinAbs(excitation_primary.get_current().value().get_processed().value().get_offset() * turnsRatioValue, max_error));

        REQUIRE_THAT(excitation_secondary.get_voltage().value().get_harmonics().value().get_amplitudes()[0], Catch::Matchers::WithinAbs(excitation_primary.get_voltage().value().get_harmonics().value().get_amplitudes()[0] /
                        turnsRatioValue, max_error));
        REQUIRE_THAT(excitation_secondary.get_voltage().value().get_harmonics().value().get_amplitudes()[1], Catch::Matchers::WithinAbs(excitation_primary.get_voltage().value().get_harmonics().value().get_amplitudes()[1] / turnsRatioValue, max_error * excitation_primary.get_voltage().value().get_harmonics().value().get_amplitudes()[1] /
                turnsRatioValue));
        REQUIRE_THAT(excitation_secondary.get_voltage().value().get_harmonics().value().get_frequencies()[0], Catch::Matchers::WithinAbs(excitation_primary.get_voltage().value().get_harmonics().value().get_frequencies()[0], max_error));
        REQUIRE_THAT(excitation_secondary.get_voltage().value().get_harmonics().value().get_frequencies()[1], Catch::Matchers::WithinAbs(excitation_primary.get_voltage().value().get_harmonics().value().get_frequencies()[1], max_error *
                        excitation_primary.get_voltage().value().get_harmonics().value().get_frequencies()[1]));

        REQUIRE_THAT(excitation_secondary.get_voltage().value().get_processed().value().get_rms().value(), Catch::Matchers::WithinAbs(excitation_primary.get_voltage().value().get_processed().value().get_rms().value() /
                        turnsRatioValue, max_error * excitation_primary.get_voltage().value().get_processed().value().get_rms().value() /
                        turnsRatioValue));
        REQUIRE_THAT(excitation_secondary.get_voltage().value().get_processed().value().get_thd().value(), Catch::Matchers::WithinAbs(excitation_primary.get_voltage().value().get_processed().value().get_thd().value(), max_error * excitation_primary.get_voltage().value().get_processed().value().get_thd().value()));
        REQUIRE_THAT(excitation_secondary.get_voltage().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(excitation_primary.get_voltage().value().get_processed().value().get_effective_frequency().value(), max_error *
                excitation_primary.get_voltage().value().get_processed().value().get_effective_frequency().value()));
        REQUIRE_THAT(excitation_secondary.get_voltage().value().get_processed().value().get_peak_to_peak().value(), Catch::Matchers::WithinAbs(excitation_primary.get_voltage().value().get_processed().value().get_peak_to_peak().value() /
                turnsRatioValue, max_error *
                excitation_primary.get_voltage().value().get_processed().value().get_peak_to_peak().value() /
                turnsRatioValue));
        REQUIRE_THAT(excitation_secondary.get_voltage().value().get_processed().value().get_offset(), Catch::Matchers::WithinAbs(excitation_primary.get_voltage().value().get_processed().value().get_offset() / turnsRatioValue, max_error));

        REQUIRE_THAT(excitation_secondary.get_voltage().value().get_harmonics().value().get_amplitudes()[0], Catch::Matchers::WithinAbs(excitation_primary.get_voltage().value().get_harmonics().value().get_amplitudes()[0] /
                        turnsRatioValue, max_error));
        REQUIRE_THAT(excitation_secondary.get_voltage().value().get_harmonics().value().get_amplitudes()[1], Catch::Matchers::WithinAbs(excitation_primary.get_voltage().value().get_harmonics().value().get_amplitudes()[1] / turnsRatioValue, max_error * excitation_primary.get_voltage().value().get_harmonics().value().get_amplitudes()[1] /
                turnsRatioValue));
        REQUIRE_THAT(excitation_secondary.get_voltage().value().get_harmonics().value().get_frequencies()[0], Catch::Matchers::WithinAbs(excitation_primary.get_voltage().value().get_harmonics().value().get_frequencies()[0], max_error *
                        excitation_primary.get_voltage().value().get_harmonics().value().get_frequencies()[0]));
        REQUIRE_THAT(excitation_secondary.get_voltage().value().get_harmonics().value().get_frequencies()[1], Catch::Matchers::WithinAbs(excitation_primary.get_voltage().value().get_harmonics().value().get_frequencies()[1], max_error *
                        excitation_primary.get_voltage().value().get_harmonics().value().get_frequencies()[1]));
    }
}

TEST_CASE("Test_One_Operating_Point_Three_Generated_Windings_Turns_Ratios_Must_Throw_Exception", "[processor][inputs][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;
    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["current"]["waveform"]["data"] = {-5, 5, -5};
    windingExcitation["current"]["waveform"]["time"] = {0, 0.0000025, 0.00001};
    windingExcitation["voltage"]["waveform"]["data"] = {1, 1, 0, 0, 1};
    windingExcitation["voltage"]["waveform"]["time"] = {0, 0.0000025, 0.0000025, 0.00001, 0.00001};
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    double turnsRatioSecondaryValue = 1;
    double turnsRatioTertiaryValue = 2;

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    json turnsRatioSecondary, turnsRatioTertiary;
    turnsRatioSecondary["nominal"] = turnsRatioSecondaryValue;
    turnsRatioTertiary["nominal"] = turnsRatioTertiaryValue;
    inputsJson["designRequirements"]["turnsRatios"] = json::array({turnsRatioSecondary, turnsRatioTertiary});

    try {
        OpenMagnetics::Inputs inputs(inputsJson);
        REQUIRE(false);
    }
    catch (const std::invalid_argument& e) {
        REQUIRE(true);
    }
}

TEST_CASE("Test_One_Operating_Point_One_Winding_Equidistant_Sinusoidal_Magnetizing_Current", "[processor][inputs][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["voltage"]["waveform"]["data"] = {0,
                                                        0.06145122063750744,
                                                        0.12290244127501487,
                                                        0.18435366191252228,
                                                        0.24580488255002975,
                                                        0.30713118165698094,
                                                        0.36841396877014293,
                                                        0.4296967558833049,
                                                        0.49097954299646696,
                                                        0.5520997341739765,
                                                        0.6130461159030302,
                                                        0.6739924976320839,
                                                        0.7349388793611376,
                                                        0.7957719836746695,
                                                        0.8562149102236266,
                                                        0.9166578367725836,
                                                        0.9771007633215409,
                                                        1.037543689870498,
                                                        1.0973397955470385,
                                                        1.1571135970555646,
                                                        1.216887398564091,
                                                        1.276661200072617,
                                                        1.3358449876276357,
                                                        1.394785828260622,
                                                        1.4537266688936081,
                                                        1.5126675095265945,
                                                        1.5711608190293265,
                                                        1.6291071460395254,
                                                        1.6870534730497244,
                                                        1.7449998000599234,
                                                        1.8027250701453021,
                                                        1.8595180566781788,
                                                        1.9163110432110557,
                                                        1.973104029743932,
                                                        2.0298970162768093,
                                                        2.085468263777627,
                                                        2.140952244204727,
                                                        2.1964362246318263,
                                                        2.251920205058926,
                                                        2.3064179538943743,
                                                        2.3604408504819907,
                                                        2.4144637470696066,
                                                        2.468486643657223,
                                                        2.5218390582265293,
                                                        2.574252797970201,
                                                        2.6266665377138723,
                                                        2.6790802774575444,
                                                        2.731216487484836,
                                                        2.78187740796716,
                                                        2.8325383284494845,
                                                        2.883199248931809,
                                                        2.933860169414133,
                                                        2.9828185803057643,
                                                        3.031587823465096,
                                                        3.0803570666244284,
                                                        3.1291263097837603,
                                                        3.176595953084751,
                                                        3.223339845815499,
                                                        3.2700837385462473,
                                                        3.3168276312769955,
                                                        3.362746026335129,
                                                        3.4073364468764473,
                                                        3.4519268674177654,
                                                        3.496517287959084,
                                                        3.540823247071316,
                                                        3.5831379761799416,
                                                        3.625452705288568,
                                                        3.6677674343971947,
                                                        3.710082163505821,
                                                        3.750324109205718,
                                                        3.790247165150427,
                                                        3.8301702210951354,
                                                        3.8700932770398437,
                                                        3.9084948307976948,
                                                        3.9459167872571586,
                                                        3.983338743716622,
                                                        4.020760700176085,
                                                        4.057271371972193,
                                                        4.092089657964925,
                                                        4.126907943957656,
                                                        4.161726229950387,
                                                        4.196297097987679,
                                                        4.228416279011962,
                                                        4.2605354600362455,
                                                        4.292654641060529,
                                                        4.324773822084812,
                                                        4.3545703852665465,
                                                        4.38390242487977,
                                                        4.413234464492994,
                                                        4.442566504106219,
                                                        4.470249709086087,
                                                        4.496714210206516,
                                                        4.523178711326947,
                                                        4.549643212447377,
                                                        4.575176689548374,
                                                        4.598701114817964,
                                                        4.6222255400875545,
                                                        4.645749965357144,
                                                        4.669099124937923,
                                                        4.689618995542171,
                                                        4.71013886614642,
                                                        4.730658736750669,
                                                        4.751178607354918,
                                                        4.7692498393956875,
                                                        4.7867089117955866,
                                                        4.8041679841954865,
                                                        4.8216270565953865,
                                                        4.8374022756616615,
                                                        4.85175269575333,
                                                        4.866103115844998,
                                                        4.880453535936666,
                                                        4.893912026714872,
                                                        4.905114460994323,
                                                        4.9163168952737735,
                                                        4.927519329553225,
                                                        4.938642296560028,
                                                        4.9466660399335725,
                                                        4.954689783307117,
                                                        4.962713526680662,
                                                        4.970737270054206,
                                                        4.976307156132361,
                                                        4.981130216077137,
                                                        4.985953276021912,
                                                        4.990776335966688,
                                                        4.993965661835207,
                                                        4.9955748186791515,
                                                        4.997183975523097,
                                                        4.998793132367041,
                                                        4.999597710789014,
                                                        4.99798855394507,
                                                        4.9963793971011246,
                                                        4.99477024025718,
                                                        4.993161083413235,
                                                        4.9883648059942995,
                                                        4.983541746049524,
                                                        4.978718686104749,
                                                        4.973895626159973,
                                                        4.9667253983674335,
                                                        4.95870165499389,
                                                        4.950677911620344,
                                                        4.942654168246801,
                                                        4.93312054669295,
                                                        4.9219181124135,
                                                        4.910715678134049,
                                                        4.899513243854598,
                                                        4.8876287459825,
                                                        4.873278325890832,
                                                        4.858927905799163,
                                                        4.844577485707495,
                                                        4.830227065615827,
                                                        4.812897520395438,
                                                        4.795438447995538,
                                                        4.777979375595638,
                                                        4.760520303195739,
                                                        4.7409186720527945,
                                                        4.720398801448546,
                                                        4.699878930844296,
                                                        4.679359060240047,
                                                        4.657512177991942,
                                                        4.633987752722351,
                                                        4.610463327452761,
                                                        4.586938902183172,
                                                        4.562875463007594,
                                                        4.536410961887164,
                                                        4.5099464607667334,
                                                        4.483481959646303,
                                                        4.457017458525873,
                                                        4.427900484299609,
                                                        4.398568444686385,
                                                        4.36923640507316,
                                                        4.339904365459936,
                                                        4.308714231572671,
                                                        4.276595050548388,
                                                        4.244475869524105,
                                                        4.212356688499822,
                                                        4.1791353729467575,
                                                        4.144317086954025,
                                                        4.109498800961293,
                                                        4.074680514968562,
                                                        4.039471678405818,
                                                        4.002049721946355,
                                                        3.96462776548689,
                                                        3.927205809027427,
                                                        3.8897838525679624,
                                                        3.85013174906749,
                                                        3.8102086931227817,
                                                        3.7702856371780733,
                                                        3.730362581233366,
                                                        3.6889247989515104,
                                                        3.6466100698428843,
                                                        3.6042953407342573,
                                                        3.5619806116256303,
                                                        3.518812498229744,
                                                        3.4742220776884256,
                                                        3.429631657147107,
                                                        3.3850412366057903,
                                                        3.3401995776423714,
                                                        3.2934556849116223,
                                                        3.246711792180875,
                                                        3.199967899450126,
                                                        3.1532240067193786,
                                                        3.1047416882040952,
                                                        3.0559724450447625,
                                                        3.0072032018854333,
                                                        2.9584339587261006,
                                                        2.9085297091729725,
                                                        2.8578687886906486,
                                                        2.807207868208325,
                                                        2.756546947726001,
                                                        2.7052871473293862,
                                                        2.652873407585714,
                                                        2.60045966784204,
                                                        2.548045928098368,
                                                        2.495498091951033,
                                                        2.44147519536342,
                                                        2.387452298775802,
                                                        2.3334294021881856,
                                                        2.2794065056005692,
                                                        2.2241782148453773,
                                                        2.1686942344182754,
                                                        2.113210253991177,
                                                        2.057726273564077,
                                                        2.001500523010371,
                                                        1.9447075364774946,
                                                        1.88791454994462,
                                                        1.8311215634117453,
                                                        1.7739729635650274,
                                                        1.7160266365548278,
                                                        1.6580803095446264,
                                                        1.6001339825344267,
                                                        1.5421379298430882,
                                                        1.483197089210103,
                                                        1.4242562485771177,
                                                        1.3653154079441308,
                                                        1.3063745673111438,
                                                        1.2467742993183535,
                                                        1.18700049780983,
                                                        1.127226696301305,
                                                        1.067452894792778,
                                                        1.0073222265960222,
                                                        0.9468793000470654,
                                                        0.8864363734981069,
                                                        0.8259934469491519,
                                                        0.7654120702256666,
                                                        0.7044656884966134,
                                                        0.6435193067675584,
                                                        0.5825729250385052,
                                                        0.5216209365530489,
                                                        0.46033814943988993,
                                                        0.39905536232672745,
                                                        0.33777257521356496,
                                                        0.27648978810040425,
                                                        0.21507927223127865,
                                                        0.15362805159377224,
                                                        0.09217683095626406,
                                                        0.03072561031875587,
                                                        -0.03072561031875054,
                                                        -0.09217683095625517,
                                                        -0.15362805159376514,
                                                        -0.21507927223126977,
                                                        -0.2764897881003954,
                                                        -0.33777257521355786,
                                                        -0.39905536232672034,
                                                        -0.4603381494398846,
                                                        -0.5216209365530453,
                                                        -0.5825729250384999,
                                                        -0.6435193067675531,
                                                        -0.7044656884966063,
                                                        -0.765412070225663,
                                                        -0.8259934469491466,
                                                        -0.8864363734981016,
                                                        -0.9468793000470601,
                                                        -1.0073222265960187,
                                                        -1.0674528947927726,
                                                        -1.1272266963012996,
                                                        -1.187000497809823,
                                                        -1.24677429931835,
                                                        -1.3063745673111367,
                                                        -1.3653154079441254,
                                                        -1.4242562485771106,
                                                        -1.4831970892100994,
                                                        -1.542137929843081,
                                                        -1.6001339825344232,
                                                        -1.6580803095446228,
                                                        -1.7160266365548225,
                                                        -1.773972963565022,
                                                        -1.8311215634117364,
                                                        -1.8879145499446164,
                                                        -1.9447075364774928,
                                                        -2.0015005230103693,
                                                        -2.057726273564077,
                                                        -2.1132102539911752,
                                                        -2.1686942344182736,
                                                        -2.224178214845372,
                                                        -2.279406505600564,
                                                        -2.333429402188182,
                                                        -2.3874522987757985,
                                                        -2.441475195363413,
                                                        -2.4954980919510312,
                                                        -2.5480459280983627,
                                                        -2.6004596678420366,
                                                        -2.6528734075857088,
                                                        -2.705287147329381,
                                                        -2.7565469477259974,
                                                        -2.8072078682083212,
                                                        -2.857868788690647,
                                                        -2.908529709172969,
                                                        -2.9584339587260953,
                                                        -3.007203201885428,
                                                        -3.055972445044759,
                                                        -3.1047416882040917,
                                                        -3.153224006719375,
                                                        -3.199967899450124,
                                                        -3.2467117921808732,
                                                        -3.2934556849116206,
                                                        -3.3401995776423696,
                                                        -3.3850412366057885,
                                                        -3.4296316571471053,
                                                        -3.474222077688424,
                                                        -3.5188124982297424,
                                                        -3.5619806116256267,
                                                        -3.6042953407342537,
                                                        -3.6466100698428807,
                                                        -3.688924798951504,
                                                        -3.7303625812333614,
                                                        -3.770285637178068,
                                                        -3.810208693122778,
                                                        -3.8501317490674865,
                                                        -3.8897838525679607,
                                                        -3.9272058090274253,
                                                        -3.96462776548689,
                                                        -4.002049721946351,
                                                        -4.039471678405816,
                                                        -4.074680514968556,
                                                        -4.109498800961289,
                                                        -4.144317086954022,
                                                        -4.179135372946753,
                                                        -4.2123566884998205,
                                                        -4.244475869524104,
                                                        -4.276595050548385,
                                                        -4.3087142315726705,
                                                        -4.339904365459933,
                                                        -4.369236405073156,
                                                        -4.398568444686381,
                                                        -4.427900484299604,
                                                        -4.45701745852587,
                                                        -4.4834819596463,
                                                        -4.509946460766729,
                                                        -4.53641096188716,
                                                        -4.56287546300759,
                                                        -4.586938902183168,
                                                        -4.610463327452758,
                                                        -4.633987752722349,
                                                        -4.657512177991937,
                                                        -4.6793590602400466,
                                                        -4.699878930844295,
                                                        -4.7203988014485425,
                                                        -4.740918672052792,
                                                        -4.760520303195736,
                                                        -4.777979375595635,
                                                        -4.795438447995536,
                                                        -4.812897520395436,
                                                        -4.830227065615827,
                                                        -4.844577485707496,
                                                        -4.858927905799163,
                                                        -4.873278325890832,
                                                        -4.887628745982499,
                                                        -4.899513243854597,
                                                        -4.910715678134047,
                                                        -4.921918112413499,
                                                        -4.933120546692948,
                                                        -4.942654168246799,
                                                        -4.950677911620344,
                                                        -4.958701654993888,
                                                        -4.9667253983674335,
                                                        -4.973895626159973,
                                                        -4.978718686104749,
                                                        -4.983541746049525,
                                                        -4.9883648059943,
                                                        -4.993161083413236,
                                                        -4.99477024025718,
                                                        -4.9963793971011246,
                                                        -4.99798855394507,
                                                        -4.999597710789014,
                                                        -4.998793132367042,
                                                        -4.997183975523097,
                                                        -4.995574818679152,
                                                        -4.993965661835208,
                                                        -4.990776335966688,
                                                        -4.985953276021913,
                                                        -4.981130216077137,
                                                        -4.9763071561323615,
                                                        -4.970737270054206,
                                                        -4.962713526680662,
                                                        -4.954689783307117,
                                                        -4.9466660399335725,
                                                        -4.938642296560029,
                                                        -4.927519329553226,
                                                        -4.916316895273774,
                                                        -4.905114460994324,
                                                        -4.893912026714872,
                                                        -4.880453535936667,
                                                        -4.866103115844997,
                                                        -4.851752695753331,
                                                        -4.837402275661661,
                                                        -4.821627056595389,
                                                        -4.804167984195489,
                                                        -4.786708911795589,
                                                        -4.769249839395689,
                                                        -4.75117860735492,
                                                        -4.730658736750671,
                                                        -4.7101388661464245,
                                                        -4.689618995542174,
                                                        -4.669099124937926,
                                                        -4.6457499653571475,
                                                        -4.622225540087557,
                                                        -4.598701114817967,
                                                        -4.5751766895483765,
                                                        -4.54964321244738,
                                                        -4.523178711326947,
                                                        -4.49671421020652,
                                                        -4.470249709086087,
                                                        -4.442566504106223,
                                                        -4.413234464493,
                                                        -4.383902424879775,
                                                        -4.35457038526655,
                                                        -4.324773822084817,
                                                        -4.2926546410605315,
                                                        -4.260535460036252,
                                                        -4.228416279011965,
                                                        -4.196297097987683,
                                                        -4.161726229950389,
                                                        -4.126907943957658,
                                                        -4.09208965796493,
                                                        -4.057271371972195,
                                                        -4.020760700176087,
                                                        -3.9833387437166223,
                                                        -3.9459167872571577,
                                                        -3.908494830797693,
                                                        -3.870093277039846,
                                                        -3.8301702210951376,
                                                        -3.7902471651504257,
                                                        -3.7503241092057173,
                                                        -3.7100821635058203,
                                                        -3.667767434397195,
                                                        -3.6254527052885734,
                                                        -3.583137976179941,
                                                        -3.5408232470713195,
                                                        -3.4965172879590867,
                                                        -3.451926867417768,
                                                        -3.407336446876453,
                                                        -3.362746026335131,
                                                        -3.316827631277004,
                                                        -3.2700837385462513,
                                                        -3.223339845815506,
                                                        -3.176595953084753,
                                                        -3.1291263097837643,
                                                        -3.0803570666244298,
                                                        -3.031587823465099,
                                                        -2.9828185803057643,
                                                        -2.9338601694141317,
                                                        -2.8831992489318097,
                                                        -2.8325383284494876,
                                                        -2.7818774079671655,
                                                        -2.73121648748484,
                                                        -2.6790802774575475,
                                                        -2.626666537713877,
                                                        -2.574252797970207,
                                                        -2.5218390582265293,
                                                        -2.4684866436572293,
                                                        -2.4144637470696075,
                                                        -2.3604408504819965,
                                                        -2.3064179538943783,
                                                        -2.25192020505893,
                                                        -2.196436224631828,
                                                        -2.1409522442047297,
                                                        -2.0854682637776314,
                                                        -2.029897016276813,
                                                        -1.9731040297439364,
                                                        -1.9163110432110635,
                                                        -1.859518056678187,
                                                        -1.802725070145307,
                                                        -1.7449998000599258,
                                                        -1.6870534730497262,
                                                        -1.62910714603953,
                                                        -1.571160819029327,
                                                        -1.512667509526601,
                                                        -1.4537266688936121,
                                                        -1.3947858282606305,
                                                        -1.3358449876276417,
                                                        -1.276661200072624,
                                                        -1.2168873985640971,
                                                        -1.1571135970555666,
                                                        -1.0973397955470396,
                                                        -1.0375436898704997,
                                                        -0.9771007633215412,
                                                        -0.9166578367725862,
                                                        -0.8562149102236276,
                                                        -0.7957719836746691,
                                                        -0.7349388793611382,
                                                        -0.6739924976320886,
                                                        -0.6130461159030389,
                                                        -0.5520997341739822,
                                                        -0.4909795429964774,
                                                        -0.4296967558833096,
                                                        -0.3684139687701524,
                                                        -0.30713118165698816,
                                                        -0.24580488255003274,
                                                        -0.18435366191252456,
                                                        -0.12290244127501992,
                                                        -0.06145122063751174,
                                                        0};
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);
    double max_error = 0.01;

    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_rms().value(), Catch::Matchers::WithinAbs(0.05627, max_error * 0.05627));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_thd().value(), Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(100000, max_error * 100000));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_peak_to_peak().value(), Catch::Matchers::WithinAbs(0.158, max_error * 0.158));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_offset(), Catch::Matchers::WithinAbs(0, max_error));

    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_harmonics().value().get_amplitudes()[0], Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_harmonics().value().get_amplitudes()[1], Catch::Matchers::WithinAbs(0.07957, max_error * 0.07957));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_harmonics().value().get_frequencies()[0], Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_harmonics().value().get_frequencies()[1], Catch::Matchers::WithinAbs(100000, max_error * 100000));
}

TEST_CASE("Test_One_Operating_Point_One_Winding_Equidistant_Triangular_Magnetizing_Current", "[processor][inputs][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["voltage"]["waveform"]["data"] = {7.5, 7.5, -2.5, -2.5, 7.5};
    windingExcitation["voltage"]["waveform"]["time"] = {0, 0.0000025, 0.0000025, 0.00001, 0.00001};
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);
    auto max_error = 0.01;

    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_rms().value(), Catch::Matchers::WithinAbs(0.0541, max_error * 0.0541));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_thd().value(), Catch::Matchers::WithinAbs(0.3765, max_error * 0.3765));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(128000, max_error * 128000));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_peak_to_peak().value(), Catch::Matchers::WithinAbs(0.1875, max_error * 0.1875));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_offset(), Catch::Matchers::WithinAbs(0, max_error));

    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_harmonics().value().get_amplitudes()[0], Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_harmonics().value().get_amplitudes()[1], Catch::Matchers::WithinAbs(0.07165, max_error * 0.07165));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_harmonics().value().get_frequencies()[0], Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_harmonics().value().get_frequencies()[1], Catch::Matchers::WithinAbs(100000, max_error * 100000));
}

TEST_CASE("Test_Quick_Operating_Point_No_Dc", "[processor][inputs][smoke-test]") {
    double dcCurrent = 0;
    double ambientTemperature = 42;
    double dutyCycle = 0.75;
    double peakToPeak = 10;
    double frequency = 100000;
    double magnetizingInductance = 10e-6;

    OpenMagnetics::Inputs inputs = OpenMagnetics::Inputs::create_quick_operating_point(
        frequency, magnetizingInductance, ambientTemperature, WaveformLabel::RECTANGULAR, peakToPeak,
        dutyCycle, dcCurrent);
    auto max_error = 0.02;

    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_rms().value(), Catch::Matchers::WithinAbs(0.541, max_error * 0.541));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_thd().value(), Catch::Matchers::WithinAbs(0.3765, max_error * 0.3765));
    // ABT #1460: exact 64-harmonic series of the piecewise-linear waveform (previously an aliased DFT128).
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(126360.81756912007, max_error * 126360.81756912007));
    // ABT #1460: exact 64-harmonic series of the piecewise-linear waveform (previously an aliased DFT128).
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_ac_effective_frequency().value(), Catch::Matchers::WithinAbs(126365.89373250987, max_error * 126365.89373250987));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_peak_to_peak().value(), Catch::Matchers::WithinAbs(1.855, max_error * 1.855));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_offset(), Catch::Matchers::WithinAbs(0, max_error));

    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_harmonics().value().get_amplitudes()[0], Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_harmonics().value().get_amplitudes()[1], Catch::Matchers::WithinAbs(0.7165, max_error * 0.7165));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_harmonics().value().get_frequencies()[0], Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_harmonics().value().get_frequencies()[1], Catch::Matchers::WithinAbs(100000, max_error * 100000));
}

TEST_CASE("Test_Quick_Operating_Point", "[processor][inputs][smoke-test]") {
    double dcCurrent = 10;
    double ambientTemperature = 42;
    double dutyCycle = 0.75;
    double peakToPeak = 10;
    double frequency = 100000;
    double magnetizingInductance = 10e-6;

    OpenMagnetics::Inputs inputs = OpenMagnetics::Inputs::create_quick_operating_point(
        frequency, magnetizingInductance, ambientTemperature, WaveformLabel::RECTANGULAR, peakToPeak,
        dutyCycle, dcCurrent);
    auto max_error = 0.01;

    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_rms().value(), Catch::Matchers::WithinAbs(10.0146, max_error * 10.0146));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_thd().value(), Catch::Matchers::WithinAbs(0.3765, max_error * 0.3765));
    // ABT #1460: exact 64-harmonic series of the piecewise-linear waveform (previously an aliased DFT128).
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(9633.2324487888418, max_error * 9633.2324487888418));
    // ABT #1460: exact 64-harmonic series of the piecewise-linear waveform (previously an aliased DFT128).
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_ac_effective_frequency().value(), Catch::Matchers::WithinAbs(126365.89373250984, max_error * 126365.89373250984));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_peak_to_peak().value(), Catch::Matchers::WithinAbs(1.855, max_error * 1.855));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_offset(), Catch::Matchers::WithinAbs(0, max_error));

    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_harmonics().value().get_amplitudes()[0], Catch::Matchers::WithinAbs(10, 10 * max_error));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_harmonics().value().get_amplitudes()[1], Catch::Matchers::WithinAbs(0.7165, max_error * 0.7165));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_harmonics().value().get_frequencies()[0], Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_harmonics().value().get_frequencies()[1], Catch::Matchers::WithinAbs(100000, max_error * 100000));
}

TEST_CASE("Test_Magnetizing_Current_Sinusoidal_Processed", "[processor][inputs][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["voltage"]["processed"]["label"] = WaveformLabel::SINUSOIDAL;
    windingExcitation["voltage"]["processed"]["offset"] = 0;
    windingExcitation["voltage"]["processed"]["peakToPeak"] = 100;
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);
    double max_error = 0.1;

    double expectedValue = 100 / 2 / (2 * M_PI * 100000 * 100e-6);

    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_peak().value(), Catch::Matchers::WithinAbs(expectedValue, max_error * expectedValue));
}

TEST_CASE("Test_Magnetizing_Current_Sinusoidal_Dc_Bias_Processed", "[processor][inputs][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["voltage"]["processed"]["label"] = WaveformLabel::SINUSOIDAL;
    windingExcitation["voltage"]["processed"]["offset"] = 0;
    windingExcitation["voltage"]["processed"]["peakToPeak"] = 100;
    windingExcitation["current"]["processed"]["label"] = WaveformLabel::SINUSOIDAL;
    windingExcitation["current"]["processed"]["offset"] = 42;
    windingExcitation["current"]["processed"]["peakToPeak"] = 1;
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);
    double max_error = 0.1;

    double expectedValue = 100 / 2 / (2 * M_PI * 100000 * 100e-6) + 42;

    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_peak().value(), Catch::Matchers::WithinAbs(expectedValue, max_error * expectedValue));
}

// Regression: pure-inductor AC excitation with V *and* I both supplied (e.g. a
// measured / pre-populated current waveform) must still produce a bipolar
// centered magnetizing current. The previous behaviour fed the excitation into
// `is_continuously_conducting_power`, which classifies V·I ≈ sin(2ωt) (V & I
// in ~90° quadrature) as discontinuous because the product spends >10% of the
// cycle near zero, then suppressed the average-subtraction of the integrated
// voltage. The integration constant of V(t) was therefore left in place,
// producing a unipolar [0, +Imax] magnetizing current and an asymmetric flux
// density whose `processed.peak` was equal to its peak-to-peak (Bpeak =
// 2 × BACpeak in the UI). See WebFrontend OP3 "HARMONICS" reproducer.
TEST_CASE("Test_Magnetizing_Current_Pure_Inductor_With_Measured_Current_Is_Bipolar", "[processor][inputs][smoke-test]") {
    json inputsJson;
    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    // Sinusoidal AC voltage, zero DC, no DC bias requested.
    json windingExcitation = json();
    windingExcitation["frequency"] = 200000;
    windingExcitation["voltage"]["processed"]["label"] = WaveformLabel::SINUSOIDAL;
    windingExcitation["voltage"]["processed"]["offset"] = 0;
    windingExcitation["voltage"]["processed"]["peakToPeak"] = 2000;
    // Pre-populated inductor current with zero DC component. This is what
    // turns on the V+I path through is_continuously_conducting_power in the
    // bug scenario.
    windingExcitation["current"]["processed"]["label"] = WaveformLabel::SINUSOIDAL;
    windingExcitation["current"]["processed"]["offset"] = 0;
    windingExcitation["current"]["processed"]["peakToPeak"] = 2;
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);
    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    REQUIRE(excitation.get_magnetizing_current().has_value());
    auto mcProcessed = excitation.get_magnetizing_current().value().get_processed().value();

    double positivePeak = mcProcessed.get_positive_peak().value();
    double negativePeak = mcProcessed.get_negative_peak().value();
    double peak = mcProcessed.get_peak().value();
    double peakToPeak = mcProcessed.get_peak_to_peak().value();
    double offset = mcProcessed.get_offset();

    // Bipolar centered: |+peak| ≈ |-peak|, midpoint ≈ 0.
    double midpoint = (positivePeak + negativePeak) / 2;
    double tolerance = 0.01 * peakToPeak;  // 1% of swing
    REQUIRE_THAT(midpoint, Catch::Matchers::WithinAbs(0.0, tolerance));
    REQUIRE_THAT(offset, Catch::Matchers::WithinAbs(0.0, tolerance));
    REQUIRE(positivePeak > 0);
    REQUIRE(negativePeak < 0);
    // peak must be the half-swing (max of |+peak|, |-peak|), not the full PP.
    REQUIRE_THAT(peak, Catch::Matchers::WithinRel(peakToPeak / 2, 0.02));
}

TEST_CASE("Test_Waveform_Coefficient_Sinusoidal", "[processor][inputs][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["voltage"]["processed"]["label"] = WaveformLabel::SINUSOIDAL;
    windingExcitation["voltage"]["processed"]["offset"] = 0;
    windingExcitation["voltage"]["processed"]["peakToPeak"] = 1000;
    windingExcitation["current"]["processed"]["label"] = WaveformLabel::SINUSOIDAL;
    windingExcitation["current"]["processed"]["offset"] = 42;
    windingExcitation["current"]["processed"]["peakToPeak"] = 1;
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);
    auto operatingPointData = inputs.get_operating_points()[0];
    double max_error = 0.1;

    double expectedValue = 4.44;

    double waveformCoefficient = OpenMagnetics::Inputs::calculate_waveform_coefficient(&operatingPointData);

    REQUIRE_THAT(waveformCoefficient, Catch::Matchers::WithinAbs(expectedValue, max_error * expectedValue));
}

TEST_CASE("Test_Waveform_Coefficient_Rectangular", "[processor][inputs][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["voltage"]["processed"]["label"] = WaveformLabel::RECTANGULAR;
    windingExcitation["voltage"]["processed"]["dutyCycle"] = 0.5;
    windingExcitation["voltage"]["processed"]["offset"] = 0;
    windingExcitation["voltage"]["processed"]["peakToPeak"] = 1000;
    windingExcitation["current"]["processed"]["label"] = WaveformLabel::SINUSOIDAL;
    windingExcitation["current"]["processed"]["offset"] = 42;
    windingExcitation["current"]["processed"]["peakToPeak"] = 1;
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);
    auto operatingPointData = inputs.get_operating_points()[0];
    double max_error = 0.1;

    double expectedValue = 4;

    double waveformCoefficient = OpenMagnetics::Inputs::calculate_waveform_coefficient(&operatingPointData);

    REQUIRE_THAT(waveformCoefficient, Catch::Matchers::WithinAbs(expectedValue, max_error * expectedValue));
}

TEST_CASE("Test_Json_Web", "[processor][inputs][smoke-test]") {

    std::string inputString = "{ \"conditions\": { \"ambientTemperature\": 42 }, \"excitationsPerWinding\": [ { \"current\": { \"harmonics\": { \"amplitudes\": [ 6.418476861114186e-16, 3.821828474081155, 1.3520347077861374, 0.4253304394765637, 1.2514181893254043e-16, 0.1536120213245558, 0.15119524913566254, 0.07875264997067084, 9.73697492775005e-17, 0.047948399955734355, 0.05513623287662622, 0.03235759328251954, 7.904827728007503e-17, 0.02339296726133173, 0.028681570420820122, 0.017771008577658842, 3.194569483808837e-16, 0.014016396066471112, 0.017807158306087294, 0.011386475028370976, 8.378597701877029e-17, 0.00947434816380318, 0.012316262725233323, 0.00804192634176095, 4.868610532679568e-17, 0.006942375707716158, 0.009173279000770505, 0.0060811614849947205, 3.6791434579594185e-17, 0.005395200034141103, 0.007217893995286599, 0.00484110135817904, 1.4304896245381993e-16, 0.00438823908696007, 0.005929248667932992, 0.004014508967799422, 8.972876846340085e-17, 0.003703646967008143, 0.005045721135369743, 0.0034434847131712696, 1.0448309804400936e-16, 0.0032247860657257607, 0.004424650744238073, 0.0030404540187283892, 8.246036444662823e-17, 0.0028849794948142492, 0.003983384372506258, 0.0027540516952759215, 1.115517517791386e-16, 0.002644278582152278, 0.0036719561660332294, 0.0025529838873121808, 7.522801894152278e-17, 0.0024780582850309095, 0.0034594522045113436, 0.002417849611852556, 1.9664496566864762e-16, 0.0023710817878943445, 0.0033268345221321993, 0.0023367953000162814, 3.474313572280248e-17, 0.002314304311502052, 0.003263064607318511, 0.002303167019445883 ], \"frequencies\": [ 0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] }, \"processed\": { \"dutyCycle\": 0.25, \"label\": \"Triangular\", \"offset\": 0, \"peakToPeak\": 10 }, \"waveform\": { \"data\": [ -5, 5, -5 ], \"time\": [ 0, 0.0000025, 0.00001 ] } }, \"frequency\": 100000, \"name\": \"My Winding Excitation\", \"voltage\": { \"harmonics\": { \"amplitudes\": [ 0, 4.502033565974343, 3.184377538608768, 1.5018839119559129, 0, 0.9025800599546525, 1.0648776307714747, 0.6462576856255092, 0, 0.5042665869042725, 0.6430565656104412, 0.4142487806223381, 0, 0.35222144081418294, 0.463801442598177, 0.30699349350336713, 0, 0.2726410924603742, 0.3654499791977766, 0.2457354325240724, 0, 0.22415467595971805, 0.3039275743318627, 0.20651574970317488, 0, 0.1918788864193534, 0.26229676419235215, 0.17958348734356613, 0, 0.16915191094100157, 0.23266758468821835, 0.16023051204701513, 0, 0.15255212918517444, 0.21087760670468925, 0.14591144183027166, 0, 0.14014832907028602, 0.19453253844864532, 0.13513635867880724, 0, 0.13077465623021037, 0.1821671775644168, 0.12698205546046792, 0, 0.12369282207450177, 0.17284496751076392, 0.12085348452158427, 0, 0.11842045801252564, 0.1659508097155576, 0.11635824698893438, 0, 0.11463807670842124, 0.1610772843508201, 0.11323684878707313, 0, 0.1121363459466216, 0.15795967496798238, 0.11132263256699867, 0, 0.11078561297814014, 0.15643843694849885, 0.1105187207328025 ], \"frequencies\": [ 0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] }, \"processed\": { \"dutyCycle\": 0.25, \"label\": \"Rectangular\", \"offset\": 0, \"peakToPeak\": 10 }, \"waveform\": { \"data\": [ 7.5, 7.5, -2.5, -2.5, 7.5 ], \"time\": [ 0, 0.0000025, 0.0000025, 0.00001, 0.00001 ] } } } ] }";
    json inputsJson = json::parse(inputString);
    OperatingPoint operatingPoint;
    from_json(inputsJson, operatingPoint);
    // Repro point: the frontend operating-point JSON must parse into a populated object.
    REQUIRE(operatingPoint.get_excitations_per_winding().size() > 0);
    CHECK(operatingPoint.get_conditions().get_ambient_temperature() == 42);
    REQUIRE(operatingPoint.get_excitations_per_winding()[0].get_current());
}

TEST_CASE("Test_Induced_Voltage", "[processor][inputs][smoke-test]") {
    double max_error = 0.1;
    double peakToPeak = 10;
    double dutyCycle = 0.25;
    double magnetizingInductance = 10e-6;
    double frequency = 100000;

    OperatingPointExcitation excitation;
    excitation.set_frequency(frequency);
    Waveform currentWaveform;
    SignalDescriptor current;
    currentWaveform.set_time(std::vector<double>{0, dutyCycle / frequency, 1.0 / frequency});
    currentWaveform.set_data(std::vector<double>{-peakToPeak / 2, peakToPeak / 2, -peakToPeak / 2});
    current.set_waveform(currentWaveform);
    excitation.set_current(current);

    double expectedValue = magnetizingInductance * frequency * peakToPeak / dutyCycle;

    auto voltage = OpenMagnetics::Inputs::calculate_induced_voltage(excitation, magnetizingInductance);

    REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(voltage.get_processed().value().get_peak().value(), max_error * expectedValue));
}

TEST_CASE("Test_Induced_Voltage_Unipolar_Triangular", "[processor][inputs][smoke-test]") {
    double max_error = 0.1;
    double peakToPeak = 10;
    double dutyCycle = 0.25;
    double magnetizingInductance = 10e-6;
    double frequency = 100000;

    OperatingPointExcitation excitation;
    excitation.set_frequency(frequency);
    Waveform currentWaveform;
    SignalDescriptor current;
    currentWaveform.set_time(std::vector<double>{0, dutyCycle / frequency, dutyCycle / frequency, 1.0 / frequency});
    currentWaveform.set_data(std::vector<double>{0, peakToPeak, 0, 0});
    current.set_waveform(currentWaveform);
    excitation.set_current(current);

    double expectedValue = magnetizingInductance * frequency * peakToPeak / dutyCycle;

    auto voltage = OpenMagnetics::Inputs::calculate_induced_voltage(excitation, magnetizingInductance);

    REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(voltage.get_processed().value().get_peak().value(), max_error * expectedValue));
}

TEST_CASE("Test_Induced_Current", "[processor][inputs][smoke-test]") {
    double max_error = 0.1;
    double peakToPeak = 53.3333;
    double dutyCycle = 0.25;
    double magnetizingInductance = 10e-6;
    double frequency = 100000;

    OperatingPointExcitation excitation;
    excitation.set_frequency(frequency);
    Waveform voltageWaveform;
    SignalDescriptor voltage;
    voltageWaveform.set_time(std::vector<double>{0, dutyCycle / frequency, dutyCycle / frequency, 1.0 / frequency, 1.0 / frequency});
    voltageWaveform.set_data(std::vector<double>{peakToPeak * (1 - dutyCycle), peakToPeak * (1 - dutyCycle), -peakToPeak * dutyCycle, -peakToPeak * dutyCycle, peakToPeak * (1 - dutyCycle)});
    voltage.set_waveform(voltageWaveform);
    excitation.set_voltage(voltage);

    double expectedValue = 5;

    auto current = OpenMagnetics::Inputs::calculate_magnetizing_current(excitation, magnetizingInductance, true, 0.0);

    REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(current.get_processed().value().get_peak().value(), max_error * expectedValue));
}

TEST_CASE("Test_Induced_Voltage_And_Magnetizing_Current_Sinusoidal", "[processor][inputs][smoke-test]") {
    Waveform waveform;
    std::vector<double> data;
    std::vector<double> time;
    double magnetizingInductance = 100e-6;
    double peakToPeak = 50;
    double frequency = 100000;
    double offset = 0;
    double max_error = 0.02;
    size_t numberPointsSampledWaveforms = 69;
    for (size_t i = 0; i < numberPointsSampledWaveforms; ++i) {
        // One full period on both axes (the angle used to step by 2*pi/69 against a time step of
        // T/68, leaving the current 2.3 A short of closing, a drift the integral then reproduced).
        double angle = i * 2 * M_PI / (numberPointsSampledWaveforms - 1);
        time.push_back(i / frequency / (numberPointsSampledWaveforms - 1));
        data.push_back((sin(angle) * peakToPeak / 2) + offset);
    }
    waveform.set_data(data);
    waveform.set_time(time);

    SignalDescriptor current;
    current.set_waveform(waveform);
    OperatingPointExcitation excitation;

    excitation.set_current(current);
    excitation.set_frequency(frequency);

    auto derivativeSignal = OpenMagnetics::Inputs::calculate_induced_voltage(excitation, magnetizingInductance);
    excitation.set_voltage(derivativeSignal);
    auto integrationSignal = OpenMagnetics::Inputs::calculate_magnetizing_current(excitation, magnetizingInductance, true, offset);

    auto calculatedProcessed = OpenMagnetics::Inputs::calculate_basic_processed_data(integrationSignal.get_waveform().value());
    REQUIRE_THAT(peakToPeak, Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * peakToPeak));
}

TEST_CASE("Test_Induced_Voltage_And_Magnetizing_Current_Bipolar_Rectangular", "[processor][inputs][smoke-test]") {
    Waveform waveform;
    std::vector<double> data;
    std::vector<double> time;
    double magnetizingInductance = 100e-6;
    double peakToPeak = 50;
    double frequency = 100000;
    double max_error = 0.02;
    waveform.set_data(std::vector<double>({0, 0, 25, 25, 0, 0, -25, -25, 0, 0}));
    waveform.set_time(std::vector<double>({0, 1.25e-6, 1.25e-6, 3.75e-6, 3.75e-6, 6.25e-6, 6.25e-6, 8.75e-6, 8.75e-6, 10e-6}));

    SignalDescriptor voltage;
    voltage.set_waveform(waveform);
    OperatingPointExcitation excitation;

    excitation.set_voltage(voltage);
    excitation.set_frequency(frequency);

    auto integrationSignal = OpenMagnetics::Inputs::calculate_magnetizing_current(excitation, magnetizingInductance, true, 0.0);
    excitation.set_current(integrationSignal);
    auto derivativeSignal = OpenMagnetics::Inputs::calculate_induced_voltage(excitation, magnetizingInductance);

    auto calculatedProcessed = OpenMagnetics::Inputs::calculate_basic_processed_data(derivativeSignal.get_waveform().value());
    REQUIRE_THAT(peakToPeak, Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * peakToPeak));
}

TEST_CASE("Test_Derivative_And_Integral_Sinusoidal", "[processor][inputs][smoke-test]") {
    Waveform waveform;
    std::vector<double> data;
    std::vector<double> time;
    double peakToPeak = 50;
    double offset = 0;
    double max_error = 0.01;
    size_t numberPointsSampledWaveforms = 69;
    for (size_t i = 0; i < numberPointsSampledWaveforms; ++i) {
        double angle = i * 2 * M_PI / numberPointsSampledWaveforms;
        time.push_back(i * 1.0 / (numberPointsSampledWaveforms - 1));
        data.push_back((sin(angle) * peakToPeak / 2) + offset);
    }
    waveform.set_data(data);
    waveform.set_time(time);

    auto integration = OpenMagnetics::Inputs::calculate_integral_waveform(waveform);
    auto derivative = OpenMagnetics::Inputs::calculate_derivative_waveform(integration);

    auto calculatedProcessed = OpenMagnetics::Inputs::calculate_basic_processed_data(derivative);
    REQUIRE_THAT(peakToPeak, Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * peakToPeak));
    REQUIRE_THAT(offset, Catch::Matchers::WithinAbs(calculatedProcessed.get_offset(), max_error));
}

TEST_CASE("Test_Derivative_And_Integral_Triangular", "[processor][inputs]") {
    SKIP("Test needs investigation");
    Waveform waveform;
    std::vector<double> data;
    std::vector<double> time;
    double peakToPeak = 10;
    double offset = 0;
    double max_error = 0.02;
    waveform.set_data(std::vector<double>({-5, 5, -5}));
    waveform.set_time(std::vector<double>({0, 0.25, 1}));

    auto derivative = OpenMagnetics::Inputs::calculate_derivative_waveform(waveform);
    auto integration = OpenMagnetics::Inputs::calculate_integral_waveform(derivative);

    auto calculatedProcessed = OpenMagnetics::Inputs::calculate_basic_processed_data(integration);
    REQUIRE_THAT(peakToPeak, Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * peakToPeak));
    REQUIRE_THAT(offset, Catch::Matchers::WithinAbs(calculatedProcessed.get_offset(), max_error));
}

TEST_CASE("Test_Derivative_And_Integral_Rectangular", "[processor][inputs][smoke-test]") {
    Waveform waveform;
    std::vector<double> data;
    std::vector<double> time;
    double peakToPeak = 10;
    double offset = 0;
    double max_error = 0.02;
    waveform.set_data(std::vector<double>({-2.5, 7.5, 7.5, -2.5, -2.5}));
    waveform.set_time(std::vector<double>({0, 0, 0.25, 0.25, 1}));

    auto integration = OpenMagnetics::Inputs::calculate_integral_waveform(waveform);
    auto derivative = OpenMagnetics::Inputs::calculate_derivative_waveform(integration);

    auto calculatedProcessed = OpenMagnetics::Inputs::calculate_basic_processed_data(derivative);
    REQUIRE_THAT(peakToPeak, Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * peakToPeak));
    REQUIRE_THAT(offset, Catch::Matchers::WithinAbs(calculatedProcessed.get_offset(), max_error));
}

TEST_CASE("Test_Derivative_And_Integral_Bipolar_Rectangular", "[processor][inputs][smoke-test]") {
    Waveform waveform;
    std::vector<double> data;
    std::vector<double> time;
    double peakToPeak = 50;
    double offset = 0;
    double max_error = 0.02;
    waveform.set_data(std::vector<double>({0, 0, 25, 25, 0, 0, -25, -25, 0, 0}));
    waveform.set_time(std::vector<double>({0, 1.25e-6, 1.25e-6, 3.75e-6, 3.75e-6, 6.25e-6, 6.25e-6, 8.75e-6, 8.75e-6, 10e-6}));

    auto integration = OpenMagnetics::Inputs::calculate_integral_waveform(waveform);

    auto derivative = OpenMagnetics::Inputs::calculate_derivative_waveform(integration);

    auto calculatedProcessed = OpenMagnetics::Inputs::calculate_basic_processed_data(derivative);
    REQUIRE_THAT(peakToPeak, Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * peakToPeak));
    REQUIRE_THAT(offset, Catch::Matchers::WithinAbs(calculatedProcessed.get_offset(), max_error));
}

TEST_CASE("Test_Reflect_Rectangular", "[processor][inputs][smoke-test]") {
    Waveform waveform;
    std::vector<double> data;
    std::vector<double> time;
    double max_error = 0.02;
    double ratio = 2;
    waveform.set_data(std::vector<double>({-2.5, 7.5, 7.5, -2.5, -2.5}));
    waveform.set_time(std::vector<double>({0, 0, 0.25, 0.25, 1}));

    SignalDescriptor signal;
    signal.set_waveform(waveform);
    
    auto reflectedSignal = OpenMagnetics::Inputs::reflect_waveform(signal, ratio);

    auto reflectedProcessed = OpenMagnetics::Inputs::calculate_basic_processed_data(reflectedSignal.get_waveform().value());
    REQUIRE_THAT(10 * ratio, Catch::Matchers::WithinAbs(reflectedProcessed.get_peak_to_peak().value(), max_error * 10 * ratio));
}

TEST_CASE("Test_Reflect_Unipolar_Rectangular_Ratio_1", "[processor][inputs][smoke-test]") {
    Waveform waveform;
    std::vector<double> data;
    std::vector<double> time;
    double max_error = 0.02;
    double ratio = 1;
    waveform.set_data(std::vector<double>({0, 60, 60, 0, 0}));
    waveform.set_time(std::vector<double>({0, 0, 0.0000025, 0.0000025, 0.00001}));

    SignalDescriptor signal;
    signal.set_waveform(waveform);
    auto processed = OpenMagnetics::Inputs::calculate_processed_data(signal.get_waveform().value(), 100000);
    
    // A unipolar pulse reflects by ampere-turn scaling, shape and timing kept (a current
    // transformer's secondary conducts during the pulse); the flyback's off-time transfer
    // belongs to the FLYBACK_* labels.
    auto reflectedSignal = OpenMagnetics::Inputs::reflect_waveform(signal, ratio, WaveformLabel::UNIPOLAR_RECTANGULAR);
    auto reflectedProcessed = OpenMagnetics::Inputs::calculate_processed_data(reflectedSignal.get_waveform().value(), 100000);

    REQUIRE_THAT(processed.get_average().value() * ratio, Catch::Matchers::WithinAbs(reflectedProcessed.get_average().value(), max_error * processed.get_average().value() * ratio));
    for (size_t i = 0; i < waveform.get_data().size(); ++i) {
        CHECK_THAT(reflectedSignal.get_waveform()->get_data()[i], Catch::Matchers::WithinAbs(waveform.get_data()[i] * ratio, 1e-12));
    }
}

TEST_CASE("Test_Reflect_Unipolar_Rectangular_Ratio_2", "[processor][inputs][smoke-test]") {
    Waveform waveform;
    std::vector<double> data;
    std::vector<double> time;
    double max_error = 0.02;
    double ratio = 2;
    waveform.set_data(std::vector<double>({0, 60, 60, 0, 0}));
    waveform.set_time(std::vector<double>({0, 0, 0.0000025, 0.0000025, 0.00001}));

    SignalDescriptor signal;
    signal.set_waveform(waveform);
    auto processed = OpenMagnetics::Inputs::calculate_processed_data(signal.get_waveform().value(), 100000);
    
    // A unipolar pulse reflects by ampere-turn scaling, shape and timing kept (a current
    // transformer's secondary conducts during the pulse); the flyback's off-time transfer
    // belongs to the FLYBACK_* labels.
    auto reflectedSignal = OpenMagnetics::Inputs::reflect_waveform(signal, ratio, WaveformLabel::UNIPOLAR_RECTANGULAR);
    auto reflectedProcessed = OpenMagnetics::Inputs::calculate_processed_data(reflectedSignal.get_waveform().value(), 100000);

    REQUIRE_THAT(processed.get_average().value() * ratio, Catch::Matchers::WithinAbs(reflectedProcessed.get_average().value(), max_error * processed.get_average().value() * ratio));
    for (size_t i = 0; i < waveform.get_data().size(); ++i) {
        CHECK_THAT(reflectedSignal.get_waveform()->get_data()[i], Catch::Matchers::WithinAbs(waveform.get_data()[i] * ratio, 1e-12));
    }
}

TEST_CASE("Test_Reflect_Unipolar_Triangular_Ratio_2", "[processor][inputs][smoke-test]") {
    Waveform waveform;
    std::vector<double> data;
    std::vector<double> time;
    double max_error = 0.02;
    double ratio = 2;
    waveform.set_data(std::vector<double>({0, 60, 0, 0}));
    waveform.set_time(std::vector<double>({0, 0.0000025, 0.0000025, 0.00001}));

    SignalDescriptor signal;
    signal.set_waveform(waveform);
    auto processed = OpenMagnetics::Inputs::calculate_processed_data(signal.get_waveform().value(), 100000);
    
    auto reflectedSignal = OpenMagnetics::Inputs::reflect_waveform(signal, ratio, WaveformLabel::UNIPOLAR_TRIANGULAR);
    auto reflectedProcessed = OpenMagnetics::Inputs::calculate_processed_data(reflectedSignal.get_waveform().value(), 100000);

    REQUIRE_THAT(processed.get_average().value() * ratio, Catch::Matchers::WithinAbs(reflectedProcessed.get_average().value(), max_error * processed.get_average().value() * ratio));
}

TEST_CASE("Test_Reflect_Flyback_Primary_Ratio_2", "[processor][inputs][smoke-test]") {
    Waveform waveform;
    std::vector<double> data;
    std::vector<double> time;
    double max_error = 0.02;
    double ratio = 2;
    waveform.set_data(std::vector<double>({0, 30, 80, 0, 0}));
    waveform.set_time(std::vector<double>({0, 0, 1.4e-6, 1.4e-6, 0.00001}));

    SignalDescriptor signal;
    signal.set_waveform(waveform);
    auto processed = OpenMagnetics::Inputs::calculate_basic_processed_data(signal.get_waveform().value());
    
    auto reflectedSignal = OpenMagnetics::Inputs::reflect_waveform(signal, ratio, WaveformLabel::FLYBACK_PRIMARY);
    auto reflectedProcessed = OpenMagnetics::Inputs::calculate_basic_processed_data(reflectedSignal.get_waveform().value());

    REQUIRE(reflectedProcessed.get_label() == WaveformLabel::FLYBACK_SECONDARY);
    REQUIRE_THAT(processed.get_peak_to_peak().value() * ratio, Catch::Matchers::WithinAbs(reflectedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value() * ratio));
}

TEST_CASE("Test_Reflect_Flyback_Secondary_Ratio_2", "[processor][inputs][smoke-test]") {
    Waveform waveform;
    std::vector<double> data;
    std::vector<double> time;
    double max_error = 0.02;
    double ratio = 2;
    waveform.set_data(std::vector<double>({0, 30, 80, 0, 0}));
    waveform.set_time(std::vector<double>({0, 0, 1.4e-6, 1.4e-6, 0.00001}));

    SignalDescriptor signal;
    signal.set_waveform(waveform);
    auto processed = OpenMagnetics::Inputs::calculate_basic_processed_data(signal.get_waveform().value());
    
    auto reflectedSignal = OpenMagnetics::Inputs::reflect_waveform(signal, ratio, WaveformLabel::FLYBACK_SECONDARY);
    auto reflectedProcessed = OpenMagnetics::Inputs::calculate_basic_processed_data(reflectedSignal.get_waveform().value());

    REQUIRE(reflectedProcessed.get_label() == WaveformLabel::FLYBACK_PRIMARY);
    REQUIRE_THAT(processed.get_peak_to_peak().value() * ratio, Catch::Matchers::WithinAbs(reflectedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value() * ratio));
}

TEST_CASE("Test_Reflect_Flyback_Primary_Web", "[processor][inputs][smoke-test]") {
    double turnRatio = 1;

    OperatingPointExcitation primaryExcitation = json::parse(R"({"current":{"harmonics":{"amplitudes":[4.558291965060757,6.529972148854835,1.789145464052662,1.942710261723393,1.1154397290543905,1.0065940903149284,0.8988582408501911,0.5788895059301848,0.7630891566656498,0.35353473619244935,0.6456068660055856,0.2676407747861538,0.5321600974207273,0.27460271311591355,0.4212435511767113,0.3071864237685926,0.31636527222157335,0.33025480179543903,0.22484933866699708,0.3342584403228846,0.1602304192861384,0.31871495949628625,0.1407347185986853,0.2864343451372125,0.16092094387433578,0.24187051785845007,0.19181915742756442,0.19085711275346715,0.21571905868290692,0.141477779475246,0.2263289212948745,0.10671871649160397,0.22229630688335306,0.10207527130577318,0.20470089894984397,0.12300075015911581,0.17626009130576112,0.14989831972782144,0.1412892563096207,0.17104729068836264,0.10662686582719509,0.1817193278103439,0.084043934810831,0.1804622345981364,0.08623921969773558,0.16772355314774284,0.10740317057400037,0.14545650387081813,0.13211333071996756,0.11727666410506915,0.15153636839094595,0.0895735293816873,0.16177942374108412,0.07382571026303117,0.16140346655602664,0.08062965411618814,0.1505782829167501,0.10267745807614781,0.1309155701859174,0.12680055354222536,0.1058030340183861,0.14569473884415238,0.08179428173758013,0.15587714970397382],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]},"processed":{"acEffectiveFrequency":621718.3722552633,"average":5.42279561360676,"dutyCycle":0.4609374999999997,"effectiveFrequency":531425.8014807898,"label":"flybackPrimary","offset":5.448133680555522,"peak":14.514900394303213,"peakToPeak":9.223090277777805,"rms":6.997011716385439,"thd":0.5671816768565314},"waveform":{"ancillaryLabel":null,"data":[0,5.448133680555522,14.671223958333327,0,0],"numberPeriods":null,"time":[0,0,0.000004609374999999998,0.000004609374999999998,0.00001]}},"frequency":100000,"magneticFieldStrength":null,"magneticFluxDensity":null,"magnetizingCurrent":null,"name":"Primary winding excitation","voltage":{"harmonics":{"amplitudes":[9.22309027777781,12.643860752984521,1.5482098237509736,3.96519064048262,1.5036230954105976,2.088180163711471,1.4309295961552153,1.194492029523459,1.332493501826718,0.6415741456288709,1.2115023045865507,0.2568357589566538,1.0718500529421924,0.024460284513222373,0.9179949822505038,0.23169804853214843,0.7547969368000894,0.3802751568040793,0.5873407469878427,0.47948718542917224,0.4207522514685872,0.5358937098386783,0.2600139243728212,0.5548945409090823,0.10978707149082016,0.5414997594159607,0.025752703324263687,0.5006963718234527,0.1430605791849741,0.4375948229458343,0.2393427194140055,0.35744706876955396,0.312647128060265,0.2655866960528321,0.3619218291128731,0.16732178127549502,0.38703829474775314,0.06780108548165185,0.3887793859782987,0.028131331656939335,0.36879243284993984,0.11608129419643629,0.3295094194096947,0.19223993817404458,0.27403749968881713,0.2535015787050763,0.2060242064300918,0.2975577510911559,0.12950268067955925,0.3229626417701206,0.048723010567191104,0.32916660187526037,0.03202370751910006,0.3165159206843967,0.10858170269138588,0.2862185126841811,0.17710037177346274,0.24027660937080772,0.23419100538841045,0.18138890635478225,0.2770642907631915,0.11282585612072807,0.3036432370169665,0.03828287613148918],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]},"processed":{"acEffectiveFrequency":593227.8773244782,"average":9.223090277777818,"dutyCycle":0.46,"effectiveFrequency":496491.712181633,"label":"custom","offset":10.004708097928473,"peak":20.009416195856947,"peakToPeak":20.009416195856947,"rms":13.584868493291284,"thd":0.49423464655132765},"waveform":{"ancillaryLabel":null,"data":[20.009416195856947,20.009416195856947,0,0],"numberPeriods":null,"time":[0,0.000004609374999999997,0.000004609374999999997,0.000009999999999999999]}}})");

    OperatingPointExcitation excitationOfThisWinding(primaryExcitation);
    auto currentSignalDescriptorProcessed = OpenMagnetics::Inputs::calculate_basic_processed_data(primaryExcitation.get_current().value().get_waveform().value());
    auto reflectedCurrentSignalDescriptor = OpenMagnetics::Inputs::reflect_waveform(primaryExcitation.get_current().value(), turnRatio, currentSignalDescriptorProcessed.get_label());
    auto processed = OpenMagnetics::Inputs::calculate_basic_processed_data(reflectedCurrentSignalDescriptor.get_waveform().value());
    REQUIRE(processed.get_label() == WaveformLabel::FLYBACK_SECONDARY);
}

TEST_CASE("Test_Reflect_Secondary_Default_Web_Rectangular_Voltage", "[processor][inputs][bug]") {
    // ABT #1670: the web's default two-winding operating point. Its primary voltage is the
    // two-level rectangle -20.5/70.5 V at duty 0.5 (with a processed block, stale, claiming
    // peakToPeak 100 and offset 0). Reflecting it 1:1 gave 0/-100 V -- a -50 V DC level across
    // a winding -- because the label guessed from the data (unipolarRectangular, the volt-
    // seconds do not balance) picked a branch that then read its numbers from the stale block.
    // A winding voltage reflects by Faraday's law: v2 = v1 / n, point by point, nothing added.
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "default_two_winding_pq2715_no_coil_found.json");
    std::ifstream file(path);
    REQUIRE(file.good());
    auto masJson = json::parse(file);
    OperatingPointExcitation primary(masJson["inputs"]["operatingPoints"][0]["excitationsPerWinding"][0]);
    auto primaryVoltage = primary.get_voltage().value();
    auto primaryData = primaryVoltage.get_waveform()->get_data();
    REQUIRE(primary.get_name().value() == "Primary winding excitation");

    for (double turnsRatio : {1.0, 2.0}) {
        auto secondary = OpenMagnetics::Inputs::calculate_reflected_secondary(primary, turnsRatio, std::string("Secondary winding excitation"));
        auto secondaryWaveform = secondary.get_voltage()->get_waveform().value();
        REQUIRE(secondaryWaveform.get_data().size() == primaryData.size());
        for (size_t i = 0; i < primaryData.size(); ++i) {
            CHECK_THAT(secondaryWaveform.get_data()[i], Catch::Matchers::WithinAbs(primaryData[i] / turnsRatio, 1e-12));
            CHECK_THAT(secondaryWaveform.get_time().value()[i], Catch::Matchers::WithinAbs(primaryVoltage.get_waveform()->get_time().value()[i], 1e-15));
        }
        // The primary's own average (+25 V, a web-side issue of its default) scales; no DC is
        // invented. Both averages go through the same sampled-waveform pipeline that fills the
        // secondary's processed block, so they are compared like for like.
        auto primarySampled = OpenMagnetics::Inputs::calculate_sampled_waveform(primaryVoltage.get_waveform().value(), primary.get_frequency());
        auto primaryVoltageReprocessed = primaryVoltage;
        primaryVoltageReprocessed.set_processed(std::nullopt);
        primaryVoltageReprocessed.set_harmonics(OpenMagnetics::Inputs::calculate_harmonics_data(primaryVoltage.get_waveform().value(), primarySampled, primary.get_frequency()));
        auto primaryAverage = OpenMagnetics::Inputs::calculate_processed_data(primaryVoltageReprocessed, primarySampled, true).get_average().value();
        CHECK(primaryAverage > 20);
        CHECK_THAT(secondary.get_voltage()->get_processed()->get_average().value(), Catch::Matchers::WithinAbs(primaryAverage / turnsRatio, 1e-9 * std::fabs(primaryAverage)));
        CHECK(secondary.get_name().value() == "Secondary winding excitation");
        // The current is a triangle: scaled by n, shape kept.
        CHECK_THAT(secondary.get_current()->get_processed()->get_peak_to_peak().value(), Catch::Matchers::WithinRel(10 * turnsRatio, 1e-6));
    }
    CHECK(!OpenMagnetics::Inputs::calculate_reflected_secondary(primary, 1.0).get_name());

    // The label-dispatched reflection (the web binding's path) must not synthesise a different
    // waveform. A rectangle straddling zero is labelled rectangular (acddb0f6) and reflects
    // by plain scaling; a unipolarRectangular that does not rest at zero has no unipolar
    // reflection.
    auto label = OpenMagnetics::Inputs::calculate_basic_processed_data(primaryVoltage.get_waveform().value()).get_label();
    REQUIRE(label == WaveformLabel::RECTANGULAR);
    auto asLabelled = OpenMagnetics::Inputs::reflect_waveform(primaryVoltage, 1.0, label);
    for (size_t i = 0; i < primaryData.size(); ++i) {
        CHECK_THAT(asLabelled.get_waveform()->get_data()[i], Catch::Matchers::WithinAbs(primaryData[i], 1e-12));
    }
    // Even forced to unipolarRectangular, the reflection is the same waveform scaled: no
    // label synthesises a shape the winding does not have.
    auto asUnipolar = OpenMagnetics::Inputs::reflect_waveform(primaryVoltage, 1.0, WaveformLabel::UNIPOLAR_RECTANGULAR);
    for (size_t i = 0; i < primaryData.size(); ++i) {
        CHECK_THAT(asUnipolar.get_waveform()->get_data()[i], Catch::Matchers::WithinAbs(primaryData[i], 1e-12));
    }
    // Shape-preserving labels reflect by plain scaling.
    auto asRectangular = OpenMagnetics::Inputs::reflect_waveform(primaryVoltage, 0.5, WaveformLabel::RECTANGULAR);
    for (size_t i = 0; i < primaryData.size(); ++i) {
        CHECK_THAT(asRectangular.get_waveform()->get_data()[i], Catch::Matchers::WithinAbs(primaryData[i] * 0.5, 1e-12));
    }
}

TEST_CASE("Test_Reflect_Secondary_Current_Transformer_Unipolar_Pulse", "[processor][inputs][bug]") {
    // A current transformer (or a forward) carries a unipolar primary current pulse, 0 A then
    // 2 A for 30% of the period. Its secondary conducts DURING the pulse: i2 = n i1, point by
    // point. The unipolar labels used to be reflected as a flyback (the current moved into
    // the off time, scaled by D/(1-D)), which is what the FLYBACK_* labels are for.
    double frequency = 100000;
    double period = 1.0 / frequency;
    double turnsRatio = 50;
    for (auto [shape, data] : std::vector<std::pair<std::string, std::vector<double>>>{
             {"rectangular", {0, 2, 2, 0, 0}}, {"triangular", {0, 2, 0, 0}}}) {
        std::vector<double> time = shape == "rectangular" ? std::vector<double>{0, 0, 0.3 * period, 0.3 * period, period}
                                                          : std::vector<double>{0, 0.3 * period, 0.3 * period, period};
        Waveform waveform;
        waveform.set_data(data);
        waveform.set_time(time);
        auto label = OpenMagnetics::Inputs::calculate_basic_processed_data(waveform).get_label();
        INFO(shape);
        CHECK(label == (shape == "rectangular" ? WaveformLabel::UNIPOLAR_RECTANGULAR : WaveformLabel::UNIPOLAR_TRIANGULAR));

        SignalDescriptor current;
        current.set_waveform(waveform);
        auto reflected = OpenMagnetics::Inputs::reflect_waveform(current, 1.0 / turnsRatio, label);
        REQUIRE(reflected.get_waveform()->get_data().size() == data.size());
        for (size_t i = 0; i < data.size(); ++i) {
            CHECK_THAT(reflected.get_waveform()->get_data()[i], Catch::Matchers::WithinAbs(data[i] / turnsRatio, 1e-12));
            CHECK_THAT(reflected.get_waveform()->get_time().value()[i], Catch::Matchers::WithinAbs(time[i], 1e-15));
        }
    }
}

TEST_CASE("Test_Reflected_Primary_Mirrors_Reflected_Secondary", "[processor][inputs][bug]") {
    // ABT #1698: the web reflected the primary in its own binding, a copy of the old logic
    // (plain current scaling, aliased harmonics). MKF owns it now: the primary from the
    // secondary is the exact inverse of calculate_reflected_secondary.
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "default_two_winding_pq2715_no_coil_found.json");
    std::ifstream file(path);
    REQUIRE(file.good());
    auto masJson = json::parse(file);
    OperatingPointExcitation primary(masJson["inputs"]["operatingPoints"][0]["excitationsPerWinding"][0]);
    auto primaryVoltageData = primary.get_voltage()->get_waveform()->get_data();
    auto primaryCurrentData = primary.get_current()->get_waveform()->get_data();

    for (double turnsRatio : {1.0, 2.0, 0.25}) {
        auto secondary = OpenMagnetics::Inputs::calculate_reflected_secondary(primary, turnsRatio, std::string("Secondary winding excitation"));
        auto back = OpenMagnetics::Inputs::calculate_reflected_primary(secondary, turnsRatio, std::string("Primary winding excitation"));
        CHECK(back.get_name().value() == "Primary winding excitation");
        auto voltage = back.get_voltage()->get_waveform()->get_data();
        auto current = back.get_current()->get_waveform()->get_data();
        REQUIRE(voltage.size() == primaryVoltageData.size());
        REQUIRE(current.size() == primaryCurrentData.size());
        for (size_t i = 0; i < voltage.size(); ++i) {
            CHECK_THAT(voltage[i], Catch::Matchers::WithinAbs(primaryVoltageData[i], 1e-9));
        }
        for (size_t i = 0; i < current.size(); ++i) {
            CHECK_THAT(current[i], Catch::Matchers::WithinAbs(primaryCurrentData[i], 1e-9));
        }
        // v1 = v2 * n: the reflected primary's processed block is computed, not copied.
        CHECK_THAT(back.get_voltage()->get_processed()->get_peak_to_peak().value(),
                   Catch::Matchers::WithinRel(secondary.get_voltage()->get_processed()->get_peak_to_peak().value() * turnsRatio, 1e-9));
        CHECK(back.get_voltage()->get_harmonics());
        CHECK(back.get_current()->get_harmonics());
    }
    CHECK(!OpenMagnetics::Inputs::calculate_reflected_primary(primary, 1.0).get_name());

    // A flyback current goes to the secondary as a flyback secondary one and comes back as a
    // flyback primary one (the label-dispatched current reflection, in both directions).
    Waveform waveform;
    waveform.set_data(std::vector<double>({0, 30, 80, 0, 0}));
    waveform.set_time(std::vector<double>({0, 0, 1.4e-6, 1.4e-6, 0.00001}));
    SignalDescriptor flybackPrimaryCurrent;
    flybackPrimaryCurrent.set_waveform(waveform);
    OperatingPointExcitation flybackPrimary(primary);
    flybackPrimary.set_current(flybackPrimaryCurrent);
    auto flybackSecondary = OpenMagnetics::Inputs::calculate_reflected_secondary(flybackPrimary, 2.0);
    CHECK(OpenMagnetics::Inputs::calculate_basic_processed_data(flybackSecondary.get_current()->get_waveform().value()).get_label() == WaveformLabel::FLYBACK_SECONDARY);
    auto flybackBack = OpenMagnetics::Inputs::calculate_reflected_primary(flybackSecondary, 2.0);
    auto flybackBackProcessed = OpenMagnetics::Inputs::calculate_basic_processed_data(flybackBack.get_current()->get_waveform().value());
    CHECK(flybackBackProcessed.get_label() == WaveformLabel::FLYBACK_PRIMARY);
    auto flybackOriginalProcessed = OpenMagnetics::Inputs::calculate_basic_processed_data(waveform);
    CHECK_THAT(flybackBackProcessed.get_peak_to_peak().value(), Catch::Matchers::WithinRel(flybackOriginalProcessed.get_peak_to_peak().value(), 0.04));

    // Missing signals throw; nothing is invented.
    OperatingPointExcitation noCurrent(primary);
    noCurrent.set_current(std::nullopt);
    REQUIRE_THROWS(OpenMagnetics::Inputs::calculate_reflected_primary(noCurrent, 1.0));
}

TEST_CASE("Test_Imported_Rectangle_Is_Custom_Not_Sinusoidal_Via_Inputs", "[processor][inputs][mas-migration][bug]") {
    // ABT #602: Inputs::calculate_basic_processed_data / try_guess_waveform_label /
    // try_guess_duty_cycle used to be a verbatim, independently-bugged twin of
    // WaveformProcessor's classifier (fixed in 6be9b794 but deliberately left
    // untouched here at the time). Every real entry point — every PyOM/WASM
    // binding — reaches THIS Inputs-facing copy, not WaveformProcessor's
    // directly, so the user-reported symptom (a 300 kHz FlyBuck switch-node
    // waveform, true duty ~20.7%, read back as label=sinusoidal, dutyCycle=0.5)
    // is only actually closed once this copy is fixed too. Reproduces the
    // reported file's shape: a two-level rectangle with finite edges and
    // plateau ripple, which compress_waveform does NOT collapse to an exact
    // 5-point analytical shape (real SPICE data never does).
    const size_t numberPoints = 512;
    const double dutyCycle = 0.207;
    const double high = 15.086;
    const double low = -3.92;
    const double period = 1.0 / 300000;
    const double ripple = 0.25;
    size_t pointsHigh = static_cast<size_t>(numberPoints * dutyCycle);

    std::vector<double> data;
    std::vector<double> time;
    for (size_t i = 0; i < numberPoints; ++i) {
        double value = (i < pointsHigh) ? high : low;
        if (i == pointsHigh || i == pointsHigh + 1) {
            value = high + (low - high) * (i - pointsHigh + 1) / 3.0;
        }
        value += ripple * sin(2 * M_PI * i * 7 / numberPoints);
        data.push_back(value);
        time.push_back(period * i / numberPoints);
    }
    Waveform waveform;
    waveform.set_data(data);
    waveform.set_time(time);

    auto processed = OpenMagnetics::Inputs::calculate_basic_processed_data(waveform);
    REQUIRE(processed.get_label() == WaveformLabel::CUSTOM);
    // Measured, not assumed: the duty must track the real high time, not be 0.5.
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(dutyCycle, 0.01));
}

TEST_CASE("Test_Flyback_Json", "[processor][inputs][smoke-test]") {
    json masJson = OpenMagneticsTesting::fixtures::get_json("mas-flyback-3-winding-lv");

    if (masJson["inputs"]["designRequirements"]["insulation"]["cti"] == "GROUP_IIIa") {
        masJson["inputs"]["designRequirements"]["insulation"]["cti"] = "GROUP_IIIA";
    }
    if (masJson["inputs"]["designRequirements"]["insulation"]["cti"] == "GROUP_IIIb") {
        masJson["inputs"]["designRequirements"]["insulation"]["cti"] = "GROUP_IIIB";
    }

    OpenMagnetics::Inputs inputs(masJson["inputs"]);
    double max_error = 0.01;

    REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
    REQUIRE_THAT(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_peak().value(), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current()->get_processed()->get_peak().value(), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current()->get_processed()->get_peak().value() * max_error));
}

TEST_CASE("Test_PFC_Json", "[processor][inputs][smoke-test]") {
    auto inventory_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "pfc_current_waveform.ndjson");
    std::ifstream ndjsonFile(inventory_path);
    std::string jsonLine;


    std::vector<double> turnsRatios;

    std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY};
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 50;
    std::getline(ndjsonFile, jsonLine);
    windingExcitation["current"]["waveform"]["time"] = json::parse(jsonLine);

    std::getline(ndjsonFile, jsonLine);
    windingExcitation["current"]["waveform"]["data"] = json::parse(jsonLine);

    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 10e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();
    OpenMagnetics::Inputs inputs(inputsJson);
    auto excitation = OpenMagnetics::Inputs::get_primary_excitation(inputs.get_mutable_operating_points()[0]);
    auto voltage = OpenMagnetics::Inputs::calculate_induced_voltage(excitation, resolve_dimensional_values(inputs.get_design_requirements().get_magnetizing_inductance()));

    excitation.set_voltage(voltage);
    inputs.get_mutable_operating_points()[0].set_excitations_per_winding({excitation});

    auto harmonics = OpenMagnetics::Inputs::get_primary_excitation(inputs.get_operating_point(0)).get_current()->get_harmonics().value();

    REQUIRE(29847U == voltage.get_waveform()->get_data().size());
    REQUIRE(22U == excitation.get_current()->get_harmonics().value().get_frequencies().size());
}

TEST_CASE("Test_Simplify_PFC_Json", "[processor][inputs]") {
    auto inventory_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "pfc_current_waveform.ndjson");
    std::ifstream ndjsonFile(inventory_path);
    std::string jsonLine;


    std::vector<double> turnsRatios;

    std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY};
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 50;
    std::getline(ndjsonFile, jsonLine);
    windingExcitation["current"]["waveform"]["time"] = json::parse(jsonLine);

    std::getline(ndjsonFile, jsonLine);
    windingExcitation["current"]["waveform"]["data"] = json::parse(jsonLine);

    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 10e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();
    double max_error = 0.03;
    {
        settings.set_inputs_trim_harmonics(false);
        OpenMagnetics::Inputs inputs(inputsJson);
        auto excitation = OpenMagnetics::Inputs::get_primary_excitation(inputs.get_mutable_operating_points()[0]);
        auto voltage = OpenMagnetics::Inputs::calculate_induced_voltage(excitation, resolve_dimensional_values(inputs.get_design_requirements().get_magnetizing_inductance()));

        excitation.set_voltage(voltage);
        inputs.get_mutable_operating_points()[0].set_excitations_per_winding({excitation});
        auto currentWaveform = excitation.get_current()->get_waveform().value();

        auto harmonics = OpenMagnetics::Inputs::get_primary_excitation(inputs.get_operating_point(0)).get_current()->get_harmonics().value();

        auto reconstructedWaveform = OpenMagnetics::Inputs::reconstruct_signal(harmonics, 50);
        auto processed = OpenMagnetics::Inputs::calculate_processed_data(harmonics, currentWaveform);
        auto reconstructedProcessed = OpenMagnetics::Inputs::calculate_processed_data(harmonics, reconstructedWaveform);


        REQUIRE_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(reconstructedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value()));
        REQUIRE_THAT(processed.get_rms().value(), Catch::Matchers::WithinAbs(reconstructedProcessed.get_rms().value(), max_error * processed.get_rms().value()));
        REQUIRE(29847U == voltage.get_waveform()->get_data().size());
        // ABT #1460: harmonics 0..N/2 of the exact piecewise-linear series (previously the FFT bins
        // 0..N/2-1), so the Nyquist bin is now included.
        REQUIRE(16385U == excitation.get_current()->get_harmonics().value().get_frequencies().size());
    }
    {
        settings.set_inputs_trim_harmonics(true);
        OpenMagnetics::Inputs inputs(inputsJson);
        auto excitation = OpenMagnetics::Inputs::get_primary_excitation(inputs.get_mutable_operating_points()[0]);
        auto voltage = OpenMagnetics::Inputs::calculate_induced_voltage(excitation, resolve_dimensional_values(inputs.get_design_requirements().get_magnetizing_inductance()));

        excitation.set_voltage(voltage);
        inputs.get_mutable_operating_points()[0].set_excitations_per_winding({excitation});
        auto currentWaveform = excitation.get_current()->get_waveform().value();

        auto harmonics = OpenMagnetics::Inputs::get_primary_excitation(inputs.get_operating_point(0)).get_current()->get_harmonics().value();

        auto reconstructedWaveform = OpenMagnetics::Inputs::reconstruct_signal(harmonics, 50);
        auto processed = OpenMagnetics::Inputs::calculate_processed_data(harmonics, currentWaveform);
        auto reconstructedProcessed = OpenMagnetics::Inputs::calculate_processed_data(harmonics, reconstructedWaveform);

        // auto outFile = outputFilePath;
        // outFile.append("Test_Simplify_PFC_Json.svg");
        // std::filesystem::remove(outFile);
        // Painter painter(outFile);
        // painter.paint_waveform(reconstructedWaveform);
        // // painter.paint_waveform(currentWaveform);
        // painter.export_svg();


        REQUIRE_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(reconstructedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value()));
        REQUIRE_THAT(processed.get_rms().value(), Catch::Matchers::WithinAbs(reconstructedProcessed.get_rms().value(), max_error * processed.get_rms().value()));
        REQUIRE(29847U == voltage.get_waveform()->get_data().size());
        REQUIRE(22U == excitation.get_current()->get_harmonics().value().get_frequencies().size());
    }
    settings.reset();
}

TEST_CASE("Test_Standardize_Waveform", "[processor][inputs][smoke-test]") {
    SignalDescriptor signalDescriptor = json::parse(R"({"harmonics":{"amplitudes":[7,3.321,10],"frequencies":[0,50,50000]},"processed":null,"waveform":null})");

    auto standardSignalDescriptor = OpenMagnetics::Inputs::standardize_waveform(signalDescriptor, 50);
    // Standardizing a harmonics-only descriptor must synthesize both a waveform and harmonics.
    REQUIRE(standardSignalDescriptor.get_waveform().has_value());
    CHECK(standardSignalDescriptor.get_waveform()->get_data().size() > 0);
    REQUIRE(standardSignalDescriptor.get_harmonics().has_value());
    CHECK(standardSignalDescriptor.get_harmonics()->get_amplitudes().size() > 0);

    #ifdef ENABLE_MATPLOTPP
        auto outFile = outputFilePath;
        outFile.append("Test_Standardize_Waveform.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_waveform(standardSignalDescriptor.get_waveform().value());
        painter.export_svg();
    #endif

}

TEST_CASE("Test_Standardize_Signal_Descriptor_Non_Multiple_Harmonics", "[processor][inputs][smoke-test]") {
    // Simulates what happens in the frontend when the user sets harmonics
    // at non-standard frequencies (e.g., 60 Hz fundamental with a 60 kHz component)
    SignalDescriptor signalDescriptor = json::parse(R"({"harmonics":{"amplitudes":[0,30.4,18.5],"frequencies":[0,60,60000]},"processed":null,"waveform":null})");

    auto standardized = OpenMagnetics::Inputs::standardize_waveform(signalDescriptor, 60);
    REQUIRE(standardized.get_waveform().has_value());
    REQUIRE(standardized.get_harmonics().has_value());

    auto processed = OpenMagnetics::Inputs::calculate_processed_data(
        standardized.get_harmonics().value(),
        standardized.get_waveform().value(), true);
    REQUIRE(processed.get_rms().has_value());
    REQUIRE(processed.get_rms().value() > 0);
}
}  // namespace

// ABT #375/#190: check_integrity resolved designRequirements.magnetizingInductance eagerly, at the
// top of the function, although exactly one branch needs it — deriving a magnetizing current for an
// excitation that carries a voltage but no current. Inputs that never reach that branch were
// rejected for a requirement nothing had asked for, with a message naming neither the field nor the
// caller: "resolve_dimensional_values: DimensionWithTolerance has neither nominal, minimum nor
// maximum set". The visible casualties were a MAS file carrying only a magnetic — from_file
// explicitly supports that, computing the inductance from the geometry and handing it to this
// constructor, but the eager resolve rejected the file before that value was ever applied — and
// Test_Flyback_Simulation.
//
// The distinction being pinned: an inputs object with NO operating points needs no inductance,
// because there is no waveform to derive anything for. One WITH an operating point does need it
// (the induced voltage of a current-driven excitation, and the magnetizing current of a
// voltage-driven one, are both proportional to L), and a missing value must still be refused there
// rather than defaulted. So this is not a relaxation — it moves the requirement to where it is real.
TEST_CASE("Test_Inputs_Magnetizing_Inductance_Resolved_Only_Where_Needed", "[processor][inputs][smoke-test]") {
    json inputsJson;
    inputsJson["operatingPoints"] = json::array();
    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"] = json::object();
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    // Nothing to process, so the absent inductance is nobody's problem.
    REQUIRE_NOTHROW(OpenMagnetics::Inputs(inputsJson));

    // Add an operating point and the inductance becomes genuinely required: it sets the induced
    // voltage across the magnetizing branch. Still absent, so this must be refused, loudly.
    json inputsWithOperatingPointJson = inputsJson;
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"]["ambientTemperature"] = 42;
    json excitation = json();
    excitation["frequency"] = 100000;
    excitation["current"]["waveform"]["data"] = {-5, 5, -5};
    excitation["current"]["waveform"]["time"] = {0, 0.0000025, 0.00001};
    operatingPoint["excitationsPerWinding"] = json::array({excitation});
    inputsWithOperatingPointJson["operatingPoints"].push_back(operatingPoint);
    REQUIRE_THROWS(OpenMagnetics::Inputs(inputsWithOperatingPointJson));

    // Supply it and the same input processes.
    json satisfiedInputsJson = inputsWithOperatingPointJson;
    satisfiedInputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    OpenMagnetics::Inputs satisfiedInputs(satisfiedInputsJson);
    REQUIRE(satisfiedInputs.get_operating_points().size() == 1);
    auto processedExcitation = satisfiedInputs.get_operating_points()[0].get_excitations_per_winding()[0];
    CHECK(processedExcitation.get_voltage());  // derived from the current and the inductance
}

TEST_CASE("Test_Reflected_Secondary_Reflects_Only_What_The_Primary_Has", "[processor][inputs]") {
    // ABT #825: with one turns ratio the design has two windings, so a lone excitation means the
    // secondary is reflected from the primary. That reflection dereferenced the primary's voltage
    // AND current unconditionally, so an operating point giving a current and no voltage — an
    // ordinary way to describe an inductor's winding — died with a bare "bad optional access"
    // naming neither the operating point, the winding, nor the field.
    json inputsJson = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array({json{{"nominal", 0.5}}});

    json excitation = json();
    excitation["frequency"] = 100000;
    excitation["current"]["processed"]["dutyCycle"] = 0.5;
    excitation["current"]["processed"]["label"] = "Triangular";
    excitation["current"]["processed"]["offset"] = 0;
    excitation["current"]["processed"]["peakToPeak"] = 10;

    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"]["ambientTemperature"] = 25;
    operatingPoint["excitationsPerWinding"] = json::array({excitation});
    inputsJson["operatingPoints"] = json::array({operatingPoint});

    // Current only: the secondary is reflected from the current alone rather than throwing.
    OpenMagnetics::Inputs inputs(inputsJson);
    auto excitations = inputs.get_operating_points()[0].get_excitations_per_winding();
    REQUIRE(excitations.size() == 2);
    CHECK(excitations[1].get_current());

    // Nothing to reflect at all: still refused, but the message must name the operating point and
    // say what is missing, never "bad optional access".
    json emptyExcitationJson = inputsJson;
    emptyExcitationJson["operatingPoints"][0]["excitationsPerWinding"][0].erase("current");
    try {
        OpenMagnetics::Inputs refused(emptyExcitationJson);
        FAIL("an excitation with neither voltage nor current must be refused");
    }
    catch (const std::exception& e) {
        std::string message(e.what());
        CHECK(message.find("Nominal") != std::string::npos);
        CHECK(message != "bad optional access");
    }
}

TEST_CASE("Test_Missing_Required_Field_Names_The_Object", "[processor][inputs]") {
    // ABT #829: a missing required field escaped as nlohmann's own
    // "[json.exception.out_of_range.403] key 'bobbin' not found", which names the key but not the
    // object — and "name" is ambiguous across the magnetic, the core, the shape, the bobbin,
    // manufacturerInfo and every winding. The field IS required (coil.json's anyOf requires it in
    // every branch); the message simply has to say whose it is.
    json coilJson = json();
    coilJson["functionalDescription"] = json::array({json{
        {"name", "primary"},
        {"numberTurns", 10},
        {"numberParallels", 1},
        {"isolationSide", "primary"},
        {"wire", "Round 0.5 - Grade 1"}}});

    try {
        auto coil = coilJson.get<OpenMagnetics::Coil>();
        FAIL("a coil without a bobbin must be refused");
    }
    catch (const std::exception& e) {
        std::string message(e.what());
        CHECK(message.find("coil") != std::string::npos);
        CHECK(message.find("bobbin") != std::string::npos);
        CHECK(message.find("out_of_range") == std::string::npos);
    }

    json windingMissingWireJson = coilJson;
    windingMissingWireJson["bobbin"] = "Dummy";
    windingMissingWireJson["functionalDescription"][0].erase("wire");
    try {
        auto coil = windingMissingWireJson.get<OpenMagnetics::Coil>();
        FAIL("a winding without a wire must be refused");
    }
    catch (const std::exception& e) {
        std::string message(e.what());
        CHECK(message.find("winding") != std::string::npos);
        CHECK(message.find("wire") != std::string::npos);
    }
}

// ABT #1330: a signal given only by processed parameters is rebuilt with create_waveform, which
// needs the duty cycle for every non-sinusoidal label and the dead time for the *_WITH_DEADTIME
// ones. When they are absent the reconstruction must refuse by name, never assume a value.
TEST_CASE("Test_Processed_Only_Signal_Defines_A_Waveform_Or_Says_What_Is_Missing", "[processor][inputs][abt-1330]") {
    ProcessedWaveform triangular;
    triangular.set_label(WaveformLabel::TRIANGULAR);
    triangular.set_peak_to_peak(1.5);
    triangular.set_offset(2.0);
    REQUIRE_THROWS_AS(OpenMagnetics::Inputs::throw_if_processed_cannot_define_waveform(triangular, "current"), InvalidInputException);
    REQUIRE_THROWS_WITH(OpenMagnetics::Inputs::throw_if_processed_cannot_define_waveform(triangular, "current"), Catch::Matchers::ContainsSubstring("dutyCycle"));
    triangular.set_duty_cycle(0.4);
    CHECK_NOTHROW(OpenMagnetics::Inputs::throw_if_processed_cannot_define_waveform(triangular, "current"));

    ProcessedWaveform withDeadTime = triangular;
    withDeadTime.set_label(WaveformLabel::RECTANGULAR_WITH_DEADTIME);
    REQUIRE_THROWS_WITH(OpenMagnetics::Inputs::throw_if_processed_cannot_define_waveform(withDeadTime, "voltage"), Catch::Matchers::ContainsSubstring("deadTime"));

    ProcessedWaveform sinusoidal;
    sinusoidal.set_label(WaveformLabel::SINUSOIDAL);
    sinusoidal.set_peak_to_peak(10);
    sinusoidal.set_offset(0);
    CHECK_NOTHROW(OpenMagnetics::Inputs::throw_if_processed_cannot_define_waveform(sinusoidal, "voltage"));
    sinusoidal.set_peak_to_peak(std::nullopt);
    REQUIRE_THROWS_WITH(OpenMagnetics::Inputs::throw_if_processed_cannot_define_waveform(sinusoidal, "voltage"), Catch::Matchers::ContainsSubstring("peakToPeak"));

    ProcessedWaveform custom = triangular;
    custom.set_label(WaveformLabel::CUSTOM);
    REQUIRE_THROWS_WITH(OpenMagnetics::Inputs::throw_if_processed_cannot_define_waveform(custom, "voltage"), Catch::Matchers::ContainsSubstring("sampled waveform"));

    // The voltage reconstruction refuses too, instead of letting create_waveform default the duty
    // cycle to 0.5: a rectangular voltage with no dutyCycle has no defined peak.
    ProcessedWaveform noDuty;
    noDuty.set_label(WaveformLabel::RECTANGULAR);
    noDuty.set_peak_to_peak(96);
    noDuty.set_offset(0);
    SignalDescriptor voltage;
    voltage.set_processed(noDuty);
    OperatingPointExcitation excitation;
    excitation.set_frequency(120000);
    excitation.set_voltage(voltage);
    OperatingPoint operatingPoint;
    operatingPoint.set_excitations_per_winding({excitation});
    OpenMagnetics::Inputs inputs;
    inputs.set_operating_points({operatingPoint});
    REQUIRE_THROWS_AS(inputs.get_maximum_voltage_peak(), InvalidInputException);
    REQUIRE_THROWS_WITH(inputs.get_maximum_voltage_peak(), Catch::Matchers::ContainsSubstring("dutyCycle"));
}

TEST_CASE("Test_Excitation_With_Proportional_Voltage_And_Current_Recalculates_Processed", "[processor][inputs][smoke-test]") {
    // The processed peak, peak to peak and offset scale with the waveform, they are not kept from the original one
    OpenMagnetics::Inputs inputs = OpenMagnetics::Inputs::create_quick_operating_point(
        100000, 100e-6, 25, WaveformLabel::RECTANGULAR, 10, 0.5, 0);
    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    auto originalVoltage = excitation.get_voltage()->get_processed().value();
    auto originalCurrent = excitation.get_current()->get_processed().value();

    auto scaled = OpenMagnetics::Inputs::get_excitation_with_proportional_voltage(excitation, 3);
    scaled = OpenMagnetics::Inputs::get_excitation_with_proportional_current(scaled, 2);
    auto scaledVoltage = scaled.get_voltage()->get_processed().value();
    auto scaledCurrent = scaled.get_current()->get_processed().value();
    REQUIRE_THAT(scaledVoltage.get_peak().value(), Catch::Matchers::WithinRel(3 * originalVoltage.get_peak().value(), 0.01));
    REQUIRE_THAT(scaledVoltage.get_peak_to_peak().value(), Catch::Matchers::WithinRel(3 * originalVoltage.get_peak_to_peak().value(), 0.01));
    REQUIRE_THAT(scaledVoltage.get_rms().value(), Catch::Matchers::WithinRel(3 * originalVoltage.get_rms().value(), 0.01));
    REQUIRE(scaledVoltage.get_label() == originalVoltage.get_label());
    REQUIRE_THAT(scaledCurrent.get_peak().value(), Catch::Matchers::WithinRel(2 * originalCurrent.get_peak().value(), 0.01));
    REQUIRE_THAT(scaledCurrent.get_peak_to_peak().value(), Catch::Matchers::WithinRel(2 * originalCurrent.get_peak_to_peak().value(), 0.01));
    REQUIRE(scaledCurrent.get_label() == originalCurrent.get_label());

    // A winding without load: no current at all
    auto unloaded = OpenMagnetics::Inputs::get_excitation_with_proportional_current(excitation, 0);
    REQUIRE_THAT(unloaded.get_current()->get_processed()->get_peak().value(), Catch::Matchers::WithinAbs(0, 1e-12));
    REQUIRE_THAT(unloaded.get_current()->get_processed()->get_rms().value(), Catch::Matchers::WithinAbs(0, 1e-12));
}

// ABT #1585: a voltage whose volt-seconds do not balance over its period has no periodic steady
// state, and its integral is a drift, not a magnetizing current. These build a rectangular winding
// voltage of +10 V for half the period and `lowVoltage` for the other half.
static Waveform abt1585_sampled_rectangular_voltage(double lowVoltage, double frequency) {
    Waveform waveform;
    double period = 1.0 / frequency;
    waveform.set_data({10, 10, lowVoltage, lowVoltage, 10});
    waveform.set_time(std::vector<double>{0, period / 2, period / 2, period, period});
    return OpenMagnetics::Inputs::calculate_sampled_waveform(waveform, frequency);
}

TEST_CASE("Test_Volt_Second_Balance_Unbalanced_Voltage_Throws", "[processor][inputs][abt-1585]") {
    double frequency = 100000;
    // Mean 2.5 V: the net 25 uV*s per period is 0.67 of the 37.5 uV*s AC swing, above the 0.25 bound.
    auto unbalanced = abt1585_sampled_rectangular_voltage(-5, frequency);
    REQUIRE_THROWS_AS(OpenMagnetics::Inputs::check_volt_second_balance(unbalanced, "test"), InvalidInputException);

    OperatingPointExcitation excitation;
    excitation.set_frequency(frequency);
    excitation.set_name("unbalanced");
    REQUIRE_THROWS_WITH(OpenMagnetics::Inputs::calculate_magnetizing_current(excitation, unbalanced, 100e-6, false),
                        Catch::Matchers::ContainsSubstring("Excitation 'unbalanced'") &&
                        Catch::Matchers::ContainsSubstring("not volt-second balanced"));
}

TEST_CASE("Test_Volt_Second_Balance_Balanced_Voltage_Integrates", "[processor][inputs][abt-1585]") {
    double frequency = 100000;
    auto balanced = abt1585_sampled_rectangular_voltage(-10, frequency);
    REQUIRE_NOTHROW(OpenMagnetics::Inputs::check_volt_second_balance(balanced, "test"));

    OperatingPointExcitation excitation;
    excitation.set_frequency(frequency);
    auto magnetizingCurrent = OpenMagnetics::Inputs::calculate_magnetizing_current(excitation, balanced, 100e-6, false);
    // 10 V over half of 10 us across 100 uH: 0.5 A peak to peak.
    REQUIRE_THAT(magnetizingCurrent.get_processed()->get_peak_to_peak().value(), Catch::Matchers::WithinRel(0.5, 0.02));
}

TEST_CASE("Test_Unbalanced_Voltage_Loads_And_Is_Reported_For_Its_Operating_Point", "[processor][inputs][abt-1855]") {
    // ABT #1855 (Alf, 2026-10-14): a voltage that is not volt-second balanced stops only its own
    // operating point's induced-current calculation (calculate_magnetizing_current of an excitation,
    // the bindings' calculate_induced_current, whose message the web shows next to the voltage).
    // Loading the design integrates its balanced part, and analyze_results reports the imbalance for
    // that operating point. Two points: a balanced +-10 V square and the +10 / -5 V one (mean 2.5 V).
    auto operatingPointJson = [](const std::string& name, double lowVoltage) {
        json operatingPoint;
        operatingPoint["name"] = name;
        operatingPoint["conditions"]["ambientTemperature"] = 25;
        json windingExcitation;
        windingExcitation["name"] = "Primary";
        windingExcitation["frequency"] = 100000;
        windingExcitation["voltage"]["waveform"]["data"] = {10, 10, lowVoltage, lowVoltage, 10};
        windingExcitation["voltage"]["waveform"]["time"] = {0, 5e-6, 5e-6, 1e-5, 1e-5};
        operatingPoint["excitationsPerWinding"] = json::array({windingExcitation});
        return operatingPoint;
    };
    json inputsJson;
    inputsJson["operatingPoints"] = json::array({operatingPointJson("Balanced", -10), operatingPointJson("Unbalanced", -5)});
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs;
    REQUIRE_NOTHROW(inputs = OpenMagnetics::Inputs(inputsJson));
    REQUIRE(inputs.get_operating_points().size() == 2);
    // The unbalanced point's magnetizing current is the integral of its balanced part: the +-7.5 V
    // square left once the 2.5 V mean is removed, 7.5 V over 5 us across 100 uH = 0.375 A peak to peak.
    auto unbalancedExcitation = inputs.get_operating_points()[1].get_excitations_per_winding()[0];
    REQUIRE(unbalancedExcitation.get_magnetizing_current());
    CHECK_THAT(unbalancedExcitation.get_magnetizing_current()->get_processed()->get_peak_to_peak().value(),
               Catch::Matchers::WithinRel(0.375, 0.02));

    CHECK(OpenMagnetics::ResultsAnalyzer::analyze_volt_second_balance(inputs.get_operating_points()[0], 0).empty());
    auto findings = OpenMagnetics::ResultsAnalyzer::analyze_volt_second_balance(inputs.get_operating_points()[1], 1);
    REQUIRE(findings.size() == 1);
    CHECK(findings[0].code == OpenMagnetics::ResultFindingCode::VOLT_SECOND_IMBALANCE);
    CHECK(findings[0].operatingPointIndex == 1);
    CHECK(findings[0].operatingPointName.value() == "Unbalanced");
    CHECK_THAT(findings[0].message, Catch::Matchers::ContainsSubstring("not volt-second balanced"));

    // The induced current asked for that voltage alone still refuses it, as the web requires.
    OperatingPointExcitation excitation = inputs.get_operating_points()[1].get_excitations_per_winding()[0];
    excitation.set_magnetizing_current(std::nullopt);
    CHECK_THROWS_WITH(OpenMagnetics::Inputs::calculate_magnetizing_current(excitation, 100e-6, false, 0.0),
                      Catch::Matchers::ContainsSubstring("not volt-second balanced"));
}

TEST_CASE("Test_Single_Winding_Inductor_With_Current_Uses_It_As_Magnetizing_Current", "[processor][inputs][abt-1585]") {
    // A line-cycle PFC-like case: the voltage is not volt-second balanced, but the winding current is
    // given, and for a single-winding inductor that current is the magnetizing current.
    json inputsJson;
    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"]["ambientTemperature"] = 25;
    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["voltage"]["waveform"]["data"] = {10, 10, -5, -5, 10};
    windingExcitation["voltage"]["waveform"]["time"] = {0, 5e-6, 5e-6, 1e-5, 1e-5};
    windingExcitation["current"]["waveform"]["data"] = {11, 13, 11};
    windingExcitation["current"]["waveform"]["time"] = {0, 5e-6, 1e-5};
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);
    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    REQUIRE(excitation.get_magnetizing_current().has_value());
    double currentPeak = excitation.get_current()->get_processed()->get_peak().value();
    double magnetizingPeak = excitation.get_magnetizing_current()->get_processed()->get_peak().value();
    REQUIRE_THAT(currentPeak, Catch::Matchers::WithinRel(13, 0.01));
    REQUIRE_THAT(magnetizingPeak, Catch::Matchers::WithinRel(currentPeak, 1e-9));
}

TEST_CASE("Test_Mas_Autocomplete_Single_Winding_Inductor_Unbalanced_Voltage_Uses_Its_Current", "[processor][inputs][abt-1585][smoke-test]") {
    // Saved web design (WebFrontend fixture etd49_wound_10uH_5T): one winding carrying a symmetric
    // +-5 A triangular current and a voltage whose volt-seconds do not balance (-20.5 / +70.5 V,
    // mean ~25 V). Since ABT #1585 the winding current of a single-winding inductor IS its
    // magnetizing current, but mas_autocomplete / inputs_autocomplete still integrated the voltage
    // and refused the whole file ("not volt-second balanced"): every such design stopped loading.
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "inductor_unbalanced_voltage_with_current.json");
    std::ifstream file(path);
    REQUIRE(file.good());
    OpenMagnetics::Mas mas(json::parse(file));
    REQUIRE(mas.get_inputs().get_design_requirements().get_turns_ratios().empty());

    auto autocompletedMas = OpenMagnetics::mas_autocomplete(mas, false, json::object());

    auto excitation = autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0];
    REQUIRE(excitation.get_magnetizing_current());
    REQUIRE(excitation.get_current());
    double magnetizingPeak = excitation.get_magnetizing_current()->get_processed()->get_peak().value();
    double currentPeak = excitation.get_current()->get_processed()->get_peak().value();
    CHECK_THAT(magnetizingPeak, Catch::Matchers::WithinRel(currentPeak, 1e-9));
    CHECK_THAT(currentPeak, Catch::Matchers::WithinRel(5.0, 1e-6));

    auto autocompletedInputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), mas.get_magnetic(), json::object());
    auto inputsExcitation = autocompletedInputs.get_operating_points()[0].get_excitations_per_winding()[0];
    CHECK_THAT(inputsExcitation.get_magnetizing_current()->get_processed()->get_peak().value(), Catch::Matchers::WithinRel(5.0, 1e-6));
}

TEST_CASE("Test_Magnetizing_Current_Removes_A_Small_Volt_Second_Imbalance", "[processor][inputs][abt-1585][smoke-test]") {
    // A measured winding voltage (oscilloscope upload) carries a probe offset, a resistive drop or a
    // capture that is not a whole number of periods. Within the tolerance the mean used to stay in
    // the integral and tilt the magnetizing current by the imbalance; a winding in periodic steady
    // state has zero net volt-seconds, so the mean is removed and the result equals the balanced one.
    auto magnetizingCurrentFor = [](double voltageOffset) {
        OperatingPointExcitation excitation;
        excitation.set_name("Primary winding excitation");
        excitation.set_frequency(100000);
        Waveform waveform;
        waveform.set_data({-50 + voltageOffset, 50 + voltageOffset, 50 + voltageOffset, -50 + voltageOffset, -50 + voltageOffset});
        waveform.set_time(std::vector<double>{0, 0, 0.000005, 0.000005, 0.00001});
        SignalDescriptor voltage;
        voltage.set_waveform(waveform);
        excitation.set_voltage(voltage);
        return OpenMagnetics::Inputs::calculate_magnetizing_current(excitation, 100e-6, false, 0.0).get_waveform()->get_data();
    };
    auto balanced = magnetizingCurrentFor(0);
    auto offset = magnetizingCurrentFor(2);  // 2 V mean: 8 % of the 25 uV*s swing, inside the 25 % limit
    REQUIRE(balanced.size() == offset.size());
    for (size_t i = 0; i < balanced.size(); ++i) {
        CHECK_THAT(offset[i], Catch::Matchers::WithinAbs(balanced[i], 1e-9));
    }
    CHECK_THROWS_AS(magnetizingCurrentFor(10), InvalidInputException);  // 40 % of the swing: still refused
}


TEST_CASE("Test_Switching_Frequency_Read_From_The_Waveform_Not_The_Cached_Harmonics", "[processor][inputs][bug]") {
    // ABT #1764: the web DMC wizard's operating point, 50 Hz line current with a 100 kHz ripple
    // sampled at 16384 points. Its processed copy (harmonics cached) and the same operating point as
    // stored (waveforms only) are the same excitation, so they have the same switching frequency.
    // get_switching_frequency returned the line frequency for the stored copy, and the magnetizing
    // inductance then saw no line swing below the switching frequency: no DC-bias roll-off, mu 300
    // instead of ~104 on the advised MPP 300 toroid, and 16.3 W of core losses instead of 1.74 W.
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "web_adviser_flows/dmc_calculate_advised_magnetics_switching_ripple.json");
    std::ifstream file(path);
    REQUIRE(file.good());
    auto fixture = json::parse(file);
    OperatingPointExcitation storedExcitation(fixture["inputs"]["operatingPoints"][0]["excitationsPerWinding"][0]);
    REQUIRE(storedExcitation.get_current()->get_waveform());
    REQUIRE_FALSE(storedExcitation.get_current()->get_harmonics());

    OpenMagnetics::Inputs inputs(fixture["inputs"]);
    auto processedExcitation = inputs.get_operating_point(0).get_excitations_per_winding()[0];
    REQUIRE(processedExcitation.get_current()->get_harmonics());

    double processedSwitchingFrequency = OpenMagnetics::Inputs::get_switching_frequency(processedExcitation);
    CHECK_THAT(processedSwitchingFrequency, Catch::Matchers::WithinRel(100000, 1e-3));
    CHECK_THAT(OpenMagnetics::Inputs::get_switching_frequency(storedExcitation), Catch::Matchers::WithinRel(processedSwitchingFrequency, 1e-9));
}

TEST_CASE("Test_Scale_Time_To_Frequency_Keeps_No_Spectrum_Of_The_Old_Time_Axis", "[processor][inputs][bug]") {
    // ABT #1764: rescaling an excitation's time axis without reprocessing kept the harmonics of the old
    // axis. The web DMC record moved to 500 Hz kept a 100 kHz ripple harmonic while its waveform's
    // ripple sat at 1 MHz. Without reprocessing the stale spectrum and processed data are dropped;
    // with it they are those of the rescaled waveform.
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "web_adviser_flows/dmc_calculate_advised_magnetics_switching_ripple.json");
    std::ifstream file(path);
    REQUIRE(file.good());
    OpenMagnetics::Inputs inputs(json::parse(file)["inputs"]);
    auto excitation = inputs.get_operating_point(0).get_excitations_per_winding()[0];
    REQUIRE(excitation.get_current()->get_harmonics());
    REQUIRE(excitation.get_voltage()->get_harmonics());

    auto unprocessed = excitation;
    OpenMagnetics::Inputs::scale_time_to_frequency(unprocessed, 500, false, false);
    CHECK_FALSE(unprocessed.get_current()->get_harmonics());
    CHECK_FALSE(unprocessed.get_current()->get_processed());
    CHECK_FALSE(unprocessed.get_voltage()->get_harmonics());
    CHECK_FALSE(unprocessed.get_voltage()->get_processed());
    REQUIRE(unprocessed.get_current()->get_waveform());

    auto reprocessed = excitation;
    OpenMagnetics::Inputs::scale_time_to_frequency(reprocessed, 500, false, true);
    auto harmonics = reprocessed.get_current()->get_harmonics().value();
    size_t strongestRipple = 2;
    for (size_t index = 2; index < harmonics.get_amplitudes().size(); ++index) {
        if (harmonics.get_amplitudes()[index] > harmonics.get_amplitudes()[strongestRipple]) {
            strongestRipple = index;
        }
    }
    CHECK_THAT(harmonics.get_frequencies()[1], Catch::Matchers::WithinRel(500, 1e-6));
    CHECK_THAT(harmonics.get_frequencies()[strongestRipple], Catch::Matchers::WithinRel(1e6, 1e-3));
}

TEST_CASE("Test_Switching_Frequency_Of_A_Line_Cycle_Record_At_Any_Frequency", "[processor][inputs][bug]") {
    // ABT #1764: a line-cycle record (a record of repeated switching cycles described at its line
    // frequency) is recognised from its waveform, not by a fixed f < 400 Hz gate. The web DMC record
    // swept to 500 Hz and 5 kHz switches at 1 MHz and 10 MHz; the gate returned 500 Hz and 5 kHz,
    // so the magnetizing inductance saw no line swing to bias the core with.
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "web_adviser_flows/dmc_calculate_advised_magnetics_switching_ripple.json");
    std::ifstream file(path);
    REQUIRE(file.good());
    auto fixture = json::parse(file);
    OperatingPointExcitation storedExcitation(fixture["inputs"]["operatingPoints"][0]["excitationsPerWinding"][0]);
    CHECK_THAT(OpenMagnetics::Inputs::get_switching_frequency(storedExcitation), Catch::Matchers::WithinRel(100000, 1e-9));
    for (double lineFrequency : {500.0, 5000.0}) {
        for (bool processSignals : {false, true}) {
            auto swept = storedExcitation;
            OpenMagnetics::Inputs::scale_time_to_frequency(swept, lineFrequency, false, processSignals);
            INFO("line frequency " << lineFrequency << " Hz, processed " << processSignals);
            CHECK_THAT(OpenMagnetics::Inputs::get_switching_frequency(swept), Catch::Matchers::WithinRel(2000 * lineFrequency, 1e-9));
        }
    }

    // One period of a converter's current is not a record of switching cycles: its frequency is
    // its switching frequency.
    auto converter = OpenMagnetics::Inputs::create_quick_operating_point(100000, 100e-6, 25, WaveformLabel::TRIANGULAR, 10, 0.5, 0);
    auto converterExcitation = converter.get_operating_point(0).get_excitations_per_winding()[0];
    CHECK_THAT(OpenMagnetics::Inputs::get_switching_frequency(converterExcitation), Catch::Matchers::WithinRel(100000, 1e-12));
}

// ABT #1467: the derived magnetizing current is Faraday on the winding voltage, never a shape
// rebuilt from the load current's label. Fixtures shared with ABT #1463 (TestWindingLosses).
namespace {
    json abt1467_inputs_without_magnetizing_current(const std::string& fileName) {
        auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), fileName);
        std::ifstream in(path);
        json masJson = json::parse(in);
        json inputsJson = masJson["inputs"];
        for (auto& operatingPoint : inputsJson["operatingPoints"]) {
            for (auto& excitation : operatingPoint["excitationsPerWinding"]) {
                excitation.erase("magnetizingCurrent");
            }
        }
        return inputsJson;
    }
}

TEST_CASE("Test_Magnetizing_Current_Dcm_Flyback_Keeps_The_Dead_Time", "[processor][inputs][abt-1467]") {
    // DCM flyback, 42:6, 100 kHz, L_m 250 uH: 235 V for 2.5 us (32 of 128 samples), -139.26 V while
    // the secondary demagnetizes (54 samples), 0 V for the last 42. The flux (and the magnetizing
    // current) rises 235 V * 2.5 us / 250 uH = 2.35 A and stays flat through the dead time. The
    // FLYBACK_PRIMARY label used to rebuild a triangle rising over D and falling over all of 1 - D:
    // no dead time. The fixture is the #1463 flyback with its reset balanced (235 * 32 / 54 V instead
    // of 137.2 V, which left 8.7 uV*s per period): any mean removed from an unbalanced voltage tilts
    // the 0 V dead time, so the dead time is flat only on a balanced one.
    OpenMagnetics::Inputs inputs(abt1467_inputs_without_magnetizing_current("Test_Magnetizing_Current_Dcm_Flyback_Balanced.json"));
    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    REQUIRE(excitation.get_magnetizing_current());
    auto magnetizingCurrent = excitation.get_magnetizing_current().value();
    auto waveform = magnetizingCurrent.get_waveform().value();
    const auto& data = waveform.get_data();
    REQUIRE(data.size() == 128);
    double maximum = *std::max_element(data.begin(), data.end());
    double minimum = *std::min_element(data.begin(), data.end());
    double faradaySwing = 235.0 * 2.5e-6 / 250e-6;
    double deadTimeLevel = data.back();
    size_t flatSamples = 0;
    for (auto value : data) {
        if (std::fabs(value - deadTimeLevel) <= 1e-9 * faradaySwing) {
            flatSamples++;
        }
    }
    INFO("i_m max " << maximum << " A, min " << minimum << " A, dead-time level " << deadTimeLevel << " A over " << flatSamples << " samples");
    CHECK(flatSamples >= 42);
    CHECK_THAT(maximum - minimum, Catch::Matchers::WithinRel(faradaySwing, 0.01));
    // DCM: the dead time sits at the flux's minimum (the balanced fixture's volt-seconds close).
    CHECK(std::fabs(deadTimeLevel - minimum) < 0.02 * faradaySwing);
    // DC anchor (ABT #907): the magnetizing current peaks where the primary current does (its
    // 128-sample peak, 2.277 A, is 3 % under the knot peak).
    CHECK_THAT(maximum, Catch::Matchers::WithinRel(faradaySwing, 0.05));
}

TEST_CASE("Test_Magnetizing_Current_Ideal_Sinusoidal_Transformer_Is_Faraday", "[processor][inputs][abt-1467]") {
    // 1:1, 20 V cos on each winding, 3 A sinusoidal currents that cancel (opposite isolation sides by
    // MKF's index convention), L_m 100 uH, 200 kHz: the magnetizing current is V / (w L_m) = 0.159 A,
    // not the 3 A load current (it used to be read as a common-mode choke and stored as a +-3 A
    // triangle).
    OpenMagnetics::Inputs inputs(abt1467_inputs_without_magnetizing_current("Test_Winding_Losses_Abt1463_Transformer.json"));
    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    REQUIRE(excitation.get_magnetizing_current());
    auto processed = excitation.get_magnetizing_current()->get_processed().value();
    double expectedAmplitude = 20.0 / (2 * std::numbers::pi * 200000 * 100e-6);
    INFO("i_m peak to peak " << processed.get_peak_to_peak().value() << " A, peak " << processed.get_peak().value() << " A");
    CHECK_THAT(processed.get_peak_to_peak().value() / 2, Catch::Matchers::WithinRel(expectedAmplitude, 0.02));
    CHECK_THAT(processed.get_peak().value(), Catch::Matchers::WithinRel(expectedAmplitude, 0.02));
}

TEST_CASE("Test_Magnetizing_Current_Of_Zero_Volt_Second_Voltage_Has_No_Dc", "[processor][inputs][abt-1467][abt-1713]") {
    // ABT #1713: 688.5 V for a quarter of the period and -229.5 V for the rest balance exactly, so
    // the flux is a zero-mean triangle. The integration constant was the trapezoid average of the
    // samples over N - 1 intervals, not the period average, leaving ~1 % of the peak as DC.
    double frequency = 100000;
    size_t numberPoints = 128;
    std::vector<double> data(numberPoints, -229.5);
    std::vector<double> time(numberPoints);
    for (size_t i = 0; i < numberPoints; ++i) {
        time[i] = i / (frequency * numberPoints);
        if (i < numberPoints / 4) {
            data[i] = 688.5;
        }
    }
    Waveform voltage;
    voltage.set_data(data);
    voltage.set_time(time);
    OperatingPointExcitation excitation;
    excitation.set_frequency(frequency);
    auto magnetizingCurrent = OpenMagnetics::Inputs::calculate_magnetizing_current(excitation, voltage, 1e-3, false, 0.0);
    double peak = magnetizingCurrent.get_processed()->get_peak().value();
    double dc = magnetizingCurrent.get_harmonics()->get_amplitudes()[0];
    double sampleMean = 0;
    auto magnetizingCurrentWaveform = magnetizingCurrent.get_waveform().value();
    for (auto value : magnetizingCurrentWaveform.get_data()) {
        sampleMean += value / numberPoints;
    }
    INFO("peak " << peak << " A, DC harmonic " << dc << " A, sample mean " << sampleMean << " A");
    REQUIRE(peak > 0);
    CHECK(std::fabs(dc) < 1e-9 * peak);
    CHECK(std::fabs(sampleMean) < 1e-9 * peak);
}

TEST_CASE("Test_Magnetizing_Current_Forward_Reset_Is_Faraday_From_Zero", "[processor][inputs][abt-1467]") {
    // Single-switch forward, 1:1, 100 kHz, judged on a core with L_m 100 uH: +100 V for 2.5 us (on),
    // -100 V for 2.5 us (reset), 0 V for 5 us (dead time). The flux rises 100 V * 2.5 us / 100 uH
    // = 2.5 A and the reset brings it back to zero, where it stays: the magnetizing current goes
    // 0 -> 2.5 A -> 0. The currents were simulated with another magnetizing inductance (200 uH:
    // the primary ramps 1.25 A over its 10 A of reflected load and resets from 1.25 A), as a
    // converter simulation judged on a candidate core always is. In the dead time no winding
    // conducts, so by Ampere the magnetizing current is zero there whatever L_m the currents had.
    // Anchoring the DC at the primary current's peak (11.25 A, load included) read an 11 A
    // magnetizing peak; at the reset's first current (1.25 A) a 1.25 A one, offset by -1.25 A
    // through the dead time; ignoring the currents (zero-mean Faraday) a 1.9 A one.
    json operatingPointJson = json::parse(R"({
        "conditions": {"ambientTemperature": 25},
        "excitationsPerWinding": [
            {"frequency": 100000, "name": "Primary",
             "current": {"waveform": {"time": [0, 2.5e-6, 2.5e-6, 5e-6, 1e-5], "data": [10, 11.25, 1.25, 0, 0]}},
             "voltage": {"waveform": {"time": [0, 2.5e-6, 2.5e-6, 5e-6, 5e-6, 1e-5], "data": [100, 100, -100, -100, 0, 0]}}},
            {"frequency": 100000, "name": "Secondary",
             "current": {"waveform": {"time": [0, 2.5e-6, 2.5e-6, 1e-5], "data": [10, 10, 0, 0]}},
             "voltage": {"waveform": {"time": [0, 2.5e-6, 2.5e-6, 5e-6, 5e-6, 1e-5], "data": [100, 100, -100, -100, 0, 0]}}}
        ]
    })");
    OperatingPoint operatingPoint(operatingPointJson);
    auto processed = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 100e-6, std::vector<double>{1});
    auto magnetizingCurrent = processed.get_excitations_per_winding()[0].get_magnetizing_current().value();
    auto waveform = magnetizingCurrent.get_waveform().value();
    const auto& data = waveform.get_data();
    double maximum = *std::max_element(data.begin(), data.end());
    double minimum = *std::min_element(data.begin(), data.end());
    INFO("i_m max " << maximum << " A, min " << minimum << " A, last sample (dead time) " << data.back() << " A");
    CHECK_THAT(maximum - minimum, Catch::Matchers::WithinRel(2.5, 0.01));
    CHECK_THAT(maximum, Catch::Matchers::WithinRel(2.5, 0.01));
    CHECK(std::fabs(minimum) < 0.01 * 2.5);
    CHECK(std::fabs(data.back()) < 0.01 * 2.5);
}

TEST_CASE("Test_Sampled_Waveform_Rejects_A_Time_Axis_That_Does_Not_Span_The_Period", "[processor][inputs][abt-1467]") {
    // A flyback primary given over its 5 us on-time only, at 100 kHz (10 us period). The sampler
    // used to take the axis length as the period and stretch the ramp over all 10 us: a 100 %
    // duty cycle and twice the RMS. It must refuse such an axis, and an axis of two periods too.
    Waveform onTimeOnly;
    onTimeOnly.set_data({0, 0, 2.35, 0});
    onTimeOnly.set_time(std::vector<double>{0, 0, 5e-6, 5e-6});
    REQUIRE_THROWS_AS(OpenMagnetics::Inputs::calculate_sampled_waveform(onTimeOnly, 100000), std::invalid_argument);

    Waveform twoPeriods;
    twoPeriods.set_data({0, 1, 0, 1, 0});
    twoPeriods.set_time(std::vector<double>{0, 5e-6, 1e-5, 1.5e-5, 2e-5});
    REQUIRE_THROWS_AS(OpenMagnetics::Inputs::calculate_sampled_waveform(twoPeriods, 100000), std::invalid_argument);

    // A frequency is required: 0 is not "the axis is the period", it throws like NaN, inf and negative.
    Waveform onePeriod;
    onePeriod.set_data({0, 1, 0});
    onePeriod.set_time(std::vector<double>{0, 5e-6, 1e-5});
    REQUIRE_THROWS_AS(OpenMagnetics::Inputs::calculate_sampled_waveform(onePeriod, 0), std::invalid_argument);
    REQUIRE_THROWS_AS(OpenMagnetics::Inputs::calculate_sampled_waveform(onePeriod, std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
    REQUIRE_THROWS_AS(OpenMagnetics::Inputs::calculate_sampled_waveform(onePeriod, std::numeric_limits<double>::infinity()), std::invalid_argument);
    REQUIRE_THROWS_AS(OpenMagnetics::Inputs::calculate_sampled_waveform(onePeriod, -100000), std::invalid_argument);
    REQUIRE_NOTHROW(OpenMagnetics::Inputs::calculate_sampled_waveform(onePeriod, 100000));

    // The same primary over the whole period, zero in the off-time, is sampled as given: the
    // second half of the samples is zero.
    Waveform fullPeriod;
    fullPeriod.set_data({0, 0, 2.35, 0, 0});
    fullPeriod.set_time(std::vector<double>{0, 0, 5e-6, 5e-6, 1e-5});
    auto sampled = OpenMagnetics::Inputs::calculate_sampled_waveform(fullPeriod, 100000);
    const auto& data = sampled.get_data();
    size_t zeroSamples = 0;
    for (auto value : data) {
        if (value == 0) {
            zeroSamples++;
        }
    }
    CHECK(zeroSamples >= data.size() / 2);
    CHECK_THAT(*std::max_element(data.begin(), data.end()), Catch::Matchers::WithinRel(2.35, 0.02));

    // A sampled waveform (open axis [0, T - dt]) spans the period too, and resampling it to the same
    // number of points gives it back: the last sample holds until T. Taking the axis length as the
    // period squeezed it by (N - 1) / N on every resampling (forward.mas.json reached 511^2 / 512^2).
    auto resampled = OpenMagnetics::Inputs::calculate_sampled_waveform(sampled, 100000, data.size());
    REQUIRE(resampled.get_data().size() == data.size());
    CHECK_THAT(resampled.get_time().value().back(), Catch::Matchers::WithinRel(sampled.get_time().value().back(), 1e-9));
    for (size_t index = 0; index < data.size(); ++index) {
        CHECK_THAT(resampled.get_data()[index], Catch::Matchers::WithinAbs(data[index], 1e-9));
    }
}

// ABT #1837: the sampled magnetizing current has a closed time axis (it reaches t0 + T). Compressed
// without the frequency it was closed a second time and spanned T + dt (113 points over 128/127 of
// the period), which the sampler then refuses: MKF rejected its own output in mas_autocomplete.
TEST_CASE("Test_Compressed_Magnetizing_Current_Spans_Exactly_One_Period", "[processor][inputs][abt-1837]") {
    double frequency = 1e6;
    double magnetizingInductance = 26e-6;
    std::vector<double> turnsRatios = {2.0};
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(frequency, magnetizingInductance, 25, WaveformLabel::SINUSOIDAL, 10, 0.5, 0, turnsRatios);
    auto operatingPoint = OpenMagnetics::Inputs::process_operating_point(inputs.get_operating_points()[0], magnetizingInductance, turnsRatios);

    for (bool compress : {false, true}) {
        CAPTURE(compress);
        auto magnetizingCurrent = OpenMagnetics::Inputs::calculate_magnetizing_current(operatingPoint, 0, magnetizingInductance, turnsRatios, compress);
        REQUIRE(magnetizingCurrent.get_waveform());
        auto time = magnetizingCurrent.get_waveform()->get_time().value();
        REQUIRE(time.size() >= 2);
        CHECK_THAT((time.back() - time.front()) * frequency, Catch::Matchers::WithinRel(1.0, 1e-9));
        // What the autocomplete does next with it.
        CHECK_NOTHROW(OpenMagnetics::Inputs::calculate_sampled_waveform(magnetizingCurrent.get_waveform().value(), frequency));
    }
}

// ABT #1878. A labelled waveform must have its label's vertices. These are the MKF flows that
// handed on a label their data did not have: the magnetizing current integrated from a
// labelled voltage copied the voltage's label (a triangle called "rectangular" or
// "bipolarRectangular"), and a current scaled to zero kept "rectangular" on 2 vertices.
TEST_CASE("Test_Magnetizing_Current_Of_A_Labelled_Voltage_Is_Not_Labelled_As_The_Voltage", "[processor][inputs][abt-1878]") {
    const double frequency = 100000;
    const double magnetizingInductance = 100e-6;
    for (auto [label, dutyCycle] : std::vector<std::pair<WaveformLabel, double>>{{WaveformLabel::RECTANGULAR, 0.3},
                                                                                 {WaveformLabel::BIPOLAR_RECTANGULAR, 0.4}}) {
        INFO(magic_enum::enum_name(label));
        auto voltageWaveform = OpenMagnetics::Inputs::create_waveform(label, 20, frequency, dutyCycle);
        REQUIRE(voltageWaveform.get_ancillary_label() == label);
        auto sampledVoltage = OpenMagnetics::Inputs::calculate_sampled_waveform(voltageWaveform, frequency);
        CHECK(!OpenMagnetics::Inputs::calculate_integral_waveform(sampledVoltage, true).get_ancillary_label());

        SignalDescriptor voltage;
        voltage.set_waveform(voltageWaveform);
        OperatingPointExcitation excitation;
        excitation.set_frequency(frequency);
        excitation.set_voltage(voltage);
        SignalDescriptor magnetizingCurrent;
        REQUIRE_NOTHROW(magnetizingCurrent = OpenMagnetics::Inputs::integrate_magnetizing_current(excitation, magnetizingInductance, true, 0.0));
        auto processed = magnetizingCurrent.get_processed().value();
        CHECK(processed.get_label() != label);
        REQUIRE(processed.get_duty_cycle());
        CHECK(processed.get_duty_cycle().value() >= 0);
        CHECK(processed.get_duty_cycle().value() <= 1);
        if (label == WaveformLabel::RECTANGULAR) {
            // A rectangle integrates to a triangle rising over the rectangle's high level.
            CHECK(processed.get_label() == WaveformLabel::TRIANGULAR);
            CHECK_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(dutyCycle, 1.0 / 128));
        }
    }
}

TEST_CASE("Test_Waveform_Scaled_By_Zero_Or_A_Negative_Factor_Keeps_No_Label", "[processor][inputs][abt-1878]") {
    const double frequency = 100000;
    auto rectangle = OpenMagnetics::Inputs::create_waveform(WaveformLabel::RECTANGULAR, 10, frequency, 0.3);
    CHECK(OpenMagnetics::Inputs::multiply_waveform(rectangle, 2).get_ancillary_label() == WaveformLabel::RECTANGULAR);
    CHECK(!OpenMagnetics::Inputs::multiply_waveform(rectangle, 0).get_ancillary_label());
    CHECK(!OpenMagnetics::Inputs::multiply_waveform(rectangle, -1).get_ancillary_label());

    OpenMagnetics::Inputs inputs = OpenMagnetics::Inputs::create_quick_operating_point(
        frequency, 100e-6, 25, WaveformLabel::RECTANGULAR, 10, 0.5, 0);
    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    OperatingPointExcitation unloaded;
    REQUIRE_NOTHROW(unloaded = OpenMagnetics::Inputs::get_excitation_with_proportional_current(excitation, 0));
    CHECK(unloaded.get_current()->get_processed()->get_label() != WaveformLabel::RECTANGULAR);
    CHECK(!unloaded.get_current()->get_waveform()->get_ancillary_label());
    OperatingPointExcitation silent;
    REQUIRE_NOTHROW(silent = OpenMagnetics::Inputs::get_excitation_with_proportional_voltage(excitation, 0));
    CHECK(!silent.get_voltage()->get_waveform()->get_ancillary_label());
}

// ABT #2022: the common-mode choke's magnetizing current is the windings' mean current. Taking the primary's own
// swing read a differential-mode line point (the two windings' currents opposite, which their harmonic amplitudes
// cannot tell from identical ones) as common-mode current: 14 A of it, 18.8 T in a 4-turn R10KZ toroid.
TEST_CASE("Common-mode choke magnetizing current is the windings' mean current", "[processor][inputs][cmc][abt-2022]") {
    auto winding_current = [](double frequency, double offset, double peakToPeak, double sign) {
        Waveform waveform;
        std::vector<double> data, time;
        size_t numberPoints = 512;
        for (size_t index = 0; index <= numberPoints; ++index) {
            double t = index / (frequency * numberPoints);
            time.push_back(t);
            data.push_back(sign * (offset + peakToPeak / 2 * std::sin(2 * std::numbers::pi * frequency * t)));
        }
        waveform.set_data(data);
        waveform.set_time(time);
        SignalDescriptor current;
        current.set_waveform(waveform);
        current.set_processed(OpenMagnetics::Inputs::calculate_processed_data(waveform, frequency, false));
        OperatingPointExcitation excitation;
        excitation.set_frequency(frequency);
        excitation.set_current(current);
        return excitation;
    };
    auto peak_to_peak = [](OperatingPoint operatingPoint) {
        return OpenMagnetics::Inputs::get_common_mode_choke_magnetizing_current(operatingPoint).get_processed()->get_peak_to_peak().value();
    };

    // Differential mode: +-14.1 A at 50 Hz, opposite in the two windings. No common-mode current drives the core.
    OperatingPoint differentialMode;
    differentialMode.set_excitations_per_winding({winding_current(50, 0, 28.28, 1), winding_current(50, 0, 28.28, -1)});
    CHECK_THAT(peak_to_peak(differentialMode), Catch::Matchers::WithinAbs(0, 1e-9));

    // Common mode on the line current: the same 1 A swing in both windings, on +10 A and -10 A of line current.
    // The line current cancels; the common-mode swing is the primary's, as before.
    OperatingPoint commonMode;
    commonMode.set_excitations_per_winding({winding_current(150000, 10, 1, 1), winding_current(150000, -10, 1, 1)});
    CHECK_THAT(peak_to_peak(commonMode), Catch::Matchers::WithinRel(1, 1e-6));
}



// ABT #2001: numbers in findings read with sensible significant figures and scaled SI units.
TEST_CASE("Test_Finding_Numbers_Use_Significant_Figures_And_SI_Prefixes", "[processor][inputs][abt-2001]") {
    CHECK(OpenMagnetics::format_quantity(1.0 / 45000, "s") == "22.2 us");
    CHECK(OpenMagnetics::format_quantity(0.000536, "V*s") == "536 uV*s");
    CHECK(OpenMagnetics::format_quantity(0.027, "Ohm") == "27 mOhm");
    CHECK(OpenMagnetics::format_quantity(999.96e-6, "s") == "1 ms");
    CHECK(OpenMagnetics::format_quantity(24.289608, "V") == "24.3 V");
    CHECK(OpenMagnetics::format_significant(26.84751, 2) == "27");
    CHECK(OpenMagnetics::format_significant(0.25) == "0.25");
    CHECK(OpenMagnetics::format_significant(2474.39) == "2474");
    CHECK(OpenMagnetics::format_significant(-0.0) == "0");
}

TEST_CASE("Test_Volt_Second_Imbalance_Message_Reads_With_Rounded_Numbers", "[processor][inputs][abt-2001]") {
    double frequency = 100000;
    // Mean 2.5 V: 24.8 uV*s against a 37.5 uV*s AC swing, 0.66 of it. calculate_volt_second_balance takes the
    // period from the sampled axis, which stops one sample short of the 10 us period (9.92 us at 128 points).
    auto unbalanced = abt1585_sampled_rectangular_voltage(-5, frequency);
    auto message = OpenMagnetics::Inputs::describe_volt_second_imbalance(unbalanced, "At operating point 0, winding 0 ('Primary')");
    REQUIRE(message);
    INFO(message.value());
    CHECK_THAT(message.value(), Catch::Matchers::StartsWith(
        "At operating point 0, winding 0 ('Primary'): the voltage is not volt-second balanced: its mean of 2.5 V drifts the "
        "flux linkage by 24.8 uV*s every 9.92 us period, 0.66 times its AC swing of 37.5 uV*s (the limit is 0.25)."));
    // No std::to_string six-decimal number is left.
    CHECK_FALSE(std::regex_search(message.value(), std::regex("[0-9]\\.[0-9]{6}")));
}

namespace {
// ABT #2002: a single-winding inductor with a large DC current, i = 900 A +- 7.5 A triangle at 45 kHz, whose
// winding voltage is L di/dt + R i: a +-VL rectangle plus the resistive drop R i(t) (24.3 V mean at 27 mOhm).
OperatingPointExcitation abt2002_resistive_inductor_excitation(double resistance, double extraOffset) {
    double frequency = 45000;
    double period = 1.0 / frequency;
    double inductance = 1.33e-6;
    double minimum = 892.5;
    double maximum = 907.5;
    double inductiveVoltage = inductance * (maximum - minimum) / (period / 2);
    Waveform current;
    current.set_time(std::vector<double>{0, period / 2, period});
    current.set_data({minimum, maximum, minimum});
    Waveform voltage;
    voltage.set_time(std::vector<double>{0, period / 2, period / 2, period});
    voltage.set_data({inductiveVoltage + resistance * minimum + extraOffset, inductiveVoltage + resistance * maximum + extraOffset,
                      -inductiveVoltage + resistance * maximum + extraOffset, -inductiveVoltage + resistance * minimum + extraOffset});
    SignalDescriptor currentSignal;
    currentSignal.set_waveform(current);
    SignalDescriptor voltageSignal;
    voltageSignal.set_waveform(voltage);
    OperatingPointExcitation excitation;
    excitation.set_frequency(frequency);
    excitation.set_name("Primary");
    excitation.set_current(currentSignal);
    excitation.set_voltage(voltageSignal);
    return excitation;
}

OperatingPoint abt2002_operating_point(const OperatingPointExcitation& excitation) {
    OperatingPoint operatingPoint;
    OperatingConditions conditions;
    conditions.set_ambient_temperature(25);
    operatingPoint.set_conditions(conditions);
    operatingPoint.set_excitations_per_winding({excitation});
    return operatingPoint;
}
}

TEST_CASE("Test_Volt_Second_Balance_Takes_Out_The_Winding_Resistive_Drop", "[processor][inputs][abt-2002]") {
    const double resistance = 0.027;
    auto excitation = abt2002_resistive_inductor_excitation(resistance, 0);
    double inductiveVoltage = 1.33e-6 * 15 / (0.5 / 45000);

    // What remains once R i(t) is taken out is the inductive rectangle, exactly, on the knots of both waveforms.
    auto magnetizing = OpenMagnetics::Inputs::remove_resistive_drop(excitation, resistance);
    auto waveform = magnetizing.get_voltage()->get_waveform().value();
    const auto& data = waveform.get_data();
    REQUIRE(data.size() == 4);
    CHECK_THAT(data[0], Catch::Matchers::WithinAbs(inductiveVoltage, 1e-9));
    CHECK_THAT(data[1], Catch::Matchers::WithinAbs(inductiveVoltage, 1e-9));
    CHECK_THAT(data[2], Catch::Matchers::WithinAbs(-inductiveVoltage, 1e-9));
    CHECK_THAT(data[3], Catch::Matchers::WithinAbs(-inductiveVoltage, 1e-9));

    auto operatingPoint = abt2002_operating_point(excitation);
    // The voltage as given: its 24.3 V mean is 27 times the AC swing, the ticket's finding.
    auto asGiven = OpenMagnetics::ResultsAnalyzer::analyze_volt_second_balance(operatingPoint, 0);
    REQUIRE(asGiven.size() == 1);
    CHECK(asGiven[0].code == OpenMagnetics::ResultFindingCode::VOLT_SECOND_IMBALANCE);
    // With the winding's resistance the mean is R * mean(i): balanced, no finding.
    CHECK(OpenMagnetics::ResultsAnalyzer::analyze_volt_second_balance(operatingPoint, 0, std::vector<double>{resistance}, 25.0).empty());

    // A genuine imbalance on top of the resistive drop is still reported, and names the drop taken out.
    auto unbalanced = abt2002_operating_point(abt2002_resistive_inductor_excitation(resistance, 5));
    auto findings = OpenMagnetics::ResultsAnalyzer::analyze_volt_second_balance(unbalanced, 0, std::vector<double>{resistance}, 25.0);
    REQUIRE(findings.size() == 1);
    INFO(findings[0].message);
    CHECK(findings[0].code == OpenMagnetics::ResultFindingCode::VOLT_SECOND_IMBALANCE);
    CHECK_THAT(findings[0].message, Catch::Matchers::ContainsSubstring("its mean of 5 V") &&
                                    Catch::Matchers::ContainsSubstring("resistive drop (27 mOhm at 25 C, 24.3 V on average)"));

    // An excitation without a current has no known drop.
    OperatingPointExcitation voltageOnly = excitation;
    voltageOnly.set_current(std::nullopt);
    CHECK_THROWS_AS(OpenMagnetics::Inputs::remove_resistive_drop(voltageOnly, resistance), InvalidInputException);
}

TEST_CASE("Test_Resistive_Drop_Of_A_Sampled_Voltage_Is_Given_Over_One_Period", "[processor][inputs][abt-2002]") {
    // A sampled voltage's axis is open, 0 .. T - dt. Built on its own knots, the corrected voltage stopped at
    // T - dt as well, a non-uniform axis that spans neither T nor T - dt once the current's knots are in it,
    // and the sampler refused it ("Waveform time axis spans 2.2e-05 s ... the excitation period is 2.22e-05 s").
    const double resistance = 0.027;
    auto excitation = abt2002_resistive_inductor_excitation(resistance, 0);
    const double frequency = excitation.get_frequency();
    const double period = 1.0 / frequency;
    const double inductiveVoltage = 1.33e-6 * 15 / (period / 2);
    const size_t numberSamples = 128;
    std::vector<double> time;
    std::vector<double> data;
    for (size_t index = 0; index < numberSamples; ++index) {
        double instant = period * index / numberSamples;
        // i(t): 892.5 A rising to 907.5 A at T/2 and back; v = +-VL + R i, the falling half from T/2 on.
        double current = instant < period / 2 ? 892.5 + 15 * instant / (period / 2) : 907.5 - 15 * (instant - period / 2) / (period / 2);
        time.push_back(instant);
        data.push_back((instant < period / 2 ? inductiveVoltage : -inductiveVoltage) + resistance * current);
    }
    Waveform sampledVoltage;
    sampledVoltage.set_time(time);
    sampledVoltage.set_data(data);
    SignalDescriptor voltageSignal;
    voltageSignal.set_waveform(sampledVoltage);
    excitation.set_voltage(voltageSignal);

    // The current's knots (0, T/2, T) are samples of the voltage: the corrected voltage is the same sampled
    // waveform, on the voltage's own open axis.
    {
        auto magnetizing = OpenMagnetics::Inputs::remove_resistive_drop(excitation, resistance);
        auto waveform = magnetizing.get_voltage()->get_waveform().value();
        const auto magnetizingTime = waveform.get_time().value();
        const auto& magnetizingData = waveform.get_data();
        CHECK(magnetizingTime == time);
        CHECK(OpenMagnetics::WaveformProcessor::check_time_axis_spans_one_period(magnetizingTime, frequency));
        REQUIRE_NOTHROW(OpenMagnetics::Inputs::calculate_sampled_waveform(waveform, frequency));
        CHECK(OpenMagnetics::Inputs::calculate_sampled_waveform(waveform, frequency).get_data().size() ==
              OpenMagnetics::Inputs::calculate_sampled_waveform(sampledVoltage, frequency).get_data().size());
        // The drop is gone on every sample: +VL on the rising half, -VL from T/2 on.
        REQUIRE(magnetizingData.size() == numberSamples);
        for (size_t index = 0; index < numberSamples; ++index) {
            INFO("sample " << index);
            CHECK_THAT(magnetizingData[index], Catch::Matchers::WithinAbs(index < numberSamples / 2 ? inductiveVoltage : -inductiveVoltage, 1e-6));
        }
    }
    // A current whose peak falls between two samples puts a knot off the voltage's grid: the corrected
    // voltage is given on the closed axis [0, T], with the knots of both.
    {
        const double peakInstant = period / 2 + period / numberSamples / 3;
        Waveform offGridCurrent;
        offGridCurrent.set_time(std::vector<double>{0, peakInstant, period});
        offGridCurrent.set_data({892.5, 907.5, 892.5});
        SignalDescriptor currentSignal;
        currentSignal.set_waveform(offGridCurrent);
        OperatingPointExcitation offGrid = excitation;
        offGrid.set_current(currentSignal);
        auto current_at = [&](double instant) {
            return instant <= peakInstant ? 892.5 + 15 * instant / peakInstant : 907.5 - 15 * (instant - peakInstant) / (period - peakInstant);
        };
        auto magnetizing = OpenMagnetics::Inputs::remove_resistive_drop(offGrid, resistance);
        auto waveform = magnetizing.get_voltage()->get_waveform().value();
        const auto magnetizingTime = waveform.get_time().value();
        const auto& magnetizingData = waveform.get_data();
        CHECK(magnetizingTime.front() == 0);
        CHECK(magnetizingTime.back() == period);
        CHECK(magnetizingTime.size() == numberSamples + 2);
        CHECK_FALSE(OpenMagnetics::WaveformProcessor::check_time_axis_spans_one_period(magnetizingTime, frequency));
        REQUIRE_NOTHROW(OpenMagnetics::Inputs::calculate_sampled_waveform(waveform, frequency));
        // Every sample of the voltage less R i at that instant; t0 + T closes back onto the first sample.
        size_t checked = 0;
        for (size_t index = 0; index < magnetizingTime.size(); ++index) {
            double samplePosition = magnetizingTime[index] / (period / numberSamples);
            if (std::fabs(samplePosition - std::round(samplePosition)) > 1e-6) {
                continue;
            }
            size_t sample = static_cast<size_t>(std::lround(samplePosition)) % numberSamples;
            INFO("sample " << samplePosition);
            CHECK_THAT(magnetizingData[index], Catch::Matchers::WithinAbs(data[sample] - resistance * current_at(magnetizingTime[index]), 1e-6));
            checked++;
        }
        CHECK(checked == numberSamples + 1);
    }
}

TEST_CASE("Test_Volt_Second_Balance_Reads_The_Wound_Coil_Resistance_At_The_Solved_Temperature", "[processor][inputs][abt-2002]") {
    auto path = std::filesystem::path{std::source_location::current().file_name()}
                    .parent_path().append("..").append("MAS").append("examples").append("02_flyback_efd25_3c95.json");
    REQUIRE(std::filesystem::exists(path));
    auto mas = OpenMagneticsTesting::mas_loader(path.string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    if (!magnetic.get_coil().get_turns_description()) {
        magnetic.get_mutable_coil().wind();
    }
    REQUIRE(magnetic.get_coil().get_turns_description());

    // No turns, no length: MKF cannot know the resistance and judges the voltage as given.
    OpenMagnetics::Coil unwound = magnetic.get_coil();
    unwound.set_turns_description(std::nullopt);
    CHECK_FALSE(OpenMagnetics::Inputs::get_winding_dc_resistances(unwound, 25));

    const double solvedTemperature = 120;
    auto resistances = OpenMagnetics::Inputs::get_winding_dc_resistances(magnetic.get_coil(), solvedTemperature);
    REQUIRE(resistances);
    auto expected = OpenMagnetics::WindingOhmicLosses::calculate_dc_resistance_per_winding(magnetic.get_coil(), solvedTemperature);
    REQUIRE(resistances->size() == expected.size());
    CHECK(resistances->at(0) == expected[0]);

    // A primary whose mean is its own drop at the solved winding temperature (the output's winding losses).
    auto excitation = abt2002_resistive_inductor_excitation(resistances->at(0), 0);
    auto operatingPoint = abt2002_operating_point(excitation);
    OpenMagnetics::Outputs output;
    WindingLossesOutput windingLosses;
    windingLosses.set_origin(ResultOrigin::SIMULATION);
    windingLosses.set_method_used("test");
    windingLosses.set_winding_losses(0);
    windingLosses.set_temperature(solvedTemperature);
    output.set_winding_losses(windingLosses);
    auto atSolved = OpenMagnetics::ResultsAnalyzer::analyze_volt_second_balance(magnetic, operatingPoint, output, 0);
    INFO((atSolved.empty() ? std::string("no finding") : atSolved[0].message));
    CHECK(atSolved.empty());
    // At the ambient (the operating point's own temperature when the output has no winding losses) the copper's
    // resistance is ~27 % lower: the remaining mean is drift and is reported.
    auto atAmbient = OpenMagnetics::ResultsAnalyzer::analyze_volt_second_balance(magnetic, operatingPoint, OpenMagnetics::Outputs(), 0);
    REQUIRE(atAmbient.size() == 1);
    CHECK_THAT(atAmbient[0].message, Catch::Matchers::ContainsSubstring("at 25 C"));
}

// ABT #1469 guards: magnetic_role, the one categorical role (gapping, the inductor loss optimiser,
// the ferrite-toroid drop, the choke handling), decided tier by tier and never throwing for missing
// descriptive metadata. The flyback is the #1467 balanced DCM flyback (42:6, 250 uH, no topology);
// the DAB is the 1.5 kW 400 V / 48 V corpus DAB (7.5:1, 490 uH, 200 kHz, no topology).
namespace {
    OpenMagnetics::Inputs abt1469_inputs(const std::string& fileName) {
        auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), fileName);
        std::ifstream in(path);
        json masJson = json::parse(in);
        return OpenMagnetics::Inputs(masJson.contains("inputs") ? masJson["inputs"] : masJson);
    }
    const std::string abt1469Flyback = "Test_Magnetizing_Current_Dcm_Flyback_Balanced.json";
    const std::string abt1469Dab = "abt1469_dab_1k5w_400v_48v_transformer.json";

    void abt1469_remove_signal(OpenMagnetics::Inputs& inputs, bool voltage, std::optional<size_t> onlyWinding = std::nullopt) {
        for (auto& operatingPoint : inputs.get_mutable_operating_points()) {
            auto& excitations = operatingPoint.get_mutable_excitations_per_winding();
            for (size_t windingIndex = 0; windingIndex < excitations.size(); ++windingIndex) {
                if (onlyWinding && windingIndex != onlyWinding.value()) {
                    continue;
                }
                if (voltage) {
                    excitations[windingIndex].set_voltage(std::nullopt);
                }
                else {
                    excitations[windingIndex].set_current(std::nullopt);
                }
            }
        }
    }
}

TEST_CASE("Test_Magnetic_Role_Declared_Tiers", "[processor][inputs][abt-1469]") {
    SECTION("One winding stores energy, whatever the topology") {
        auto inputs = abt1469_inputs(abt1469Flyback);
        inputs.get_mutable_design_requirements().set_turns_ratios({});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::nullopt);
        inputs.get_mutable_design_requirements().set_topology(MAS::Topology::LLC_RESONANT_CONVERTER);
        for (auto& operatingPoint : inputs.get_mutable_operating_points()) {
            operatingPoint.get_mutable_excitations_per_winding().resize(1);
        }
        REQUIRE(OpenMagnetics::is_single_winding_design(inputs));
        auto decision = OpenMagnetics::magnetic_role(inputs);
        CHECK(decision.role == OpenMagnetics::MagneticRole::STORING);
        CHECK(decision.tier == OpenMagnetics::MagneticRoleTier::SINGLE_WINDING);
        // A choke topology names the component: a one-winding differential-mode choke stays one, for
        // every site that reads the declared choke (the DMC core pre-filter, its repetitions).
        inputs.get_mutable_design_requirements().set_topology(MAS::Topology::DIFFERENTIAL_MODE_CHOKE);
        decision = OpenMagnetics::magnetic_role(inputs);
        CHECK(decision.role == OpenMagnetics::MagneticRole::DIFFERENTIAL_MODE_CHOKE);
        CHECK(decision.tier == OpenMagnetics::MagneticRoleTier::TOPOLOGY);
        // Its subApplication alone does not: one winding stores energy, as before the role existed.
        inputs.get_mutable_design_requirements().set_topology(std::nullopt);
        inputs.get_mutable_design_requirements().set_sub_application("differentialModeNoiseFiltering");
        CHECK(OpenMagnetics::magnetic_role(inputs).tier == OpenMagnetics::MagneticRoleTier::SINGLE_WINDING);
    }
    SECTION("The topology table decides a multi-winding design") {
        auto inputs = abt1469_inputs(abt1469Flyback);
        REQUIRE_FALSE(OpenMagnetics::is_single_winding_design(inputs));
        std::vector<std::pair<MAS::Topology, OpenMagnetics::MagneticRole>> expected{
            {MAS::Topology::DUAL_ACTIVE_BRIDGE_CONVERTER, OpenMagnetics::MagneticRole::TRANSFORMER},
            {MAS::Topology::FLYBACK_CONVERTER, OpenMagnetics::MagneticRole::STORING},
            {MAS::Topology::COMMON_MODE_CHOKE, OpenMagnetics::MagneticRole::COMMON_MODE_CHOKE},
            {MAS::Topology::DIFFERENTIAL_MODE_CHOKE, OpenMagnetics::MagneticRole::DIFFERENTIAL_MODE_CHOKE}};
        for (auto [topology, role] : expected) {
            INFO(magic_enum::enum_name(topology));
            inputs.get_mutable_design_requirements().set_topology(topology);
            auto decision = OpenMagnetics::magnetic_role(inputs);
            CHECK(decision.role == role);
            CHECK(decision.tier == OpenMagnetics::MagneticRoleTier::TOPOLOGY);
        }
        // A topology whose model returns either of two magnetics (the Weinberg converter) does not
        // decide: the next tier does.
        inputs.get_mutable_design_requirements().set_topology(MAS::Topology::WEINBERG_CONVERTER);
        CHECK_FALSE(OpenMagnetics::declared_magnetic_role(inputs));
        CHECK(OpenMagnetics::magnetic_role(inputs).tier == OpenMagnetics::MagneticRoleTier::STORAGE_RATIO);
    }
    SECTION("The subApplication decides the chokes and power filtering, not transforming") {
        auto inputs = abt1469_inputs(abt1469Dab);
        REQUIRE_FALSE(inputs.get_design_requirements().get_topology());
        std::vector<std::pair<std::string, OpenMagnetics::MagneticRole>> expected{
            {"commonModeNoiseFiltering", OpenMagnetics::MagneticRole::COMMON_MODE_CHOKE},
            {"differentialModeNoiseFiltering", OpenMagnetics::MagneticRole::DIFFERENTIAL_MODE_CHOKE},
            {"powerFiltering", OpenMagnetics::MagneticRole::STORING}};
        for (auto [subApplication, role] : expected) {
            INFO(subApplication);
            inputs.get_mutable_design_requirements().set_sub_application(subApplication);
            auto decision = OpenMagnetics::magnetic_role(inputs);
            CHECK(decision.role == role);
            CHECK(decision.tier == OpenMagnetics::MagneticRoleTier::APPLICATION);
        }
        // The subApplication-declared differential-mode choke is what the magnetizing model is told.
        inputs.get_mutable_design_requirements().set_sub_application("differentialModeNoiseFiltering");
        CHECK(OpenMagnetics::magnetizing_model_topology(inputs) == MAS::Topology::DIFFERENTIAL_MODE_CHOKE);
        // A flyback isolates and transforms too: neither says transformer.
        for (std::string subApplication : {"transforming", "isolation"}) {
            INFO(subApplication);
            auto flyback = abt1469_inputs(abt1469Flyback);
            flyback.get_mutable_design_requirements().set_sub_application(subApplication);
            CHECK_FALSE(OpenMagnetics::declared_magnetic_role(flyback));
            CHECK(OpenMagnetics::magnetic_role(flyback).role == OpenMagnetics::MagneticRole::STORING);
            CHECK_FALSE(OpenMagnetics::magnetizing_model_topology(flyback));
        }
    }
    SECTION("Every winding on one isolation side stores energy") {
        auto inputs = abt1469_inputs(abt1469Dab);
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY, IsolationSide::PRIMARY});
        auto decision = OpenMagnetics::magnetic_role(inputs);
        CHECK(decision.role == OpenMagnetics::MagneticRole::STORING);
        CHECK(decision.tier == OpenMagnetics::MagneticRoleTier::ISOLATION_SIDES);
    }
}

TEST_CASE("Test_Magnetic_Role_From_The_Waveforms", "[processor][inputs][abt-1469]") {
    SECTION("Storage ratio: the DAB transfers, the flyback stores") {
        auto dab = abt1469_inputs(abt1469Dab);
        REQUIRE_FALSE(OpenMagnetics::declared_magnetic_role(dab));
        auto dabDecision = OpenMagnetics::magnetic_role(dab);
        INFO("DAB storage ratio " << dabDecision.storageRatio.value_or(-1));
        CHECK(dabDecision.tier == OpenMagnetics::MagneticRoleTier::STORAGE_RATIO);
        CHECK(dabDecision.role == OpenMagnetics::MagneticRole::TRANSFORMER);
        REQUIRE(dabDecision.storageRatio);
        CHECK(dabDecision.storageRatio.value() < OpenMagnetics::magneticRoleStoringThreshold);

        auto flyback = abt1469_inputs(abt1469Flyback);
        REQUIRE_FALSE(OpenMagnetics::declared_magnetic_role(flyback));
        auto flybackDecision = OpenMagnetics::magnetic_role(flyback);
        INFO("flyback storage ratio " << flybackDecision.storageRatio.value_or(-1));
        CHECK(flybackDecision.tier == OpenMagnetics::MagneticRoleTier::STORAGE_RATIO);
        CHECK(flybackDecision.role == OpenMagnetics::MagneticRole::STORING);
        REQUIRE(flybackDecision.storageRatio);
        CHECK(flybackDecision.storageRatio.value() >= OpenMagnetics::magneticRoleStoringThreshold);
    }
    SECTION("Net magnetomotive force, from the currents alone") {
        auto dab = abt1469_inputs(abt1469Dab);
        abt1469_remove_signal(dab, /*voltage=*/true);
        auto dabDecision = OpenMagnetics::magnetic_role(dab);
        INFO("DAB net MMF ratio " << dabDecision.netMagnetomotiveForceRatio.value_or(-1));
        CHECK(dabDecision.tier == OpenMagnetics::MagneticRoleTier::NET_MAGNETOMOTIVE_FORCE);
        CHECK(dabDecision.role == OpenMagnetics::MagneticRole::TRANSFORMER);

        auto flyback = abt1469_inputs(abt1469Flyback);
        abt1469_remove_signal(flyback, /*voltage=*/true);
        auto flybackDecision = OpenMagnetics::magnetic_role(flyback);
        INFO("flyback net MMF ratio " << flybackDecision.netMagnetomotiveForceRatio.value_or(-1));
        CHECK(flybackDecision.tier == OpenMagnetics::MagneticRoleTier::NET_MAGNETOMOTIVE_FORCE);
        CHECK(flybackDecision.role == OpenMagnetics::MagneticRole::STORING);
    }
    SECTION("Voltages only, or one winding's current missing: a transformer, without throwing") {
        auto voltagesOnly = abt1469_inputs(abt1469Flyback);
        abt1469_remove_signal(voltagesOnly, /*voltage=*/false);
        OpenMagnetics::MagneticRoleDecision decision;
        REQUIRE_NOTHROW(decision = OpenMagnetics::magnetic_role(voltagesOnly));
        CHECK(decision.tier == OpenMagnetics::MagneticRoleTier::NO_CURRENTS);
        CHECK(decision.role == OpenMagnetics::MagneticRole::TRANSFORMER);

        auto partialCurrents = abt1469_inputs(abt1469Flyback);
        abt1469_remove_signal(partialCurrents, /*voltage=*/false, 1);
        REQUIRE_NOTHROW(decision = OpenMagnetics::magnetic_role(partialCurrents));
        CHECK(decision.tier == OpenMagnetics::MagneticRoleTier::NO_CURRENTS);
        CHECK(decision.role == OpenMagnetics::MagneticRole::TRANSFORMER);
    }
}

TEST_CASE("Test_Max_Volt_Seconds_Is_The_Peak_Flux_Linkage", "[processor][inputs][abt-2070]") {
    // ABT #2070: calculate_max_volt_seconds feeds B_peak = lambda / (N * A_e) and the transformer
    // turn seed N = lambda / (B_max * A_e), so it must be the PEAK flux linkage. It used to return
    // max|int_0^t V dt|, which for a sine starting at 0 V is the whole excursion 2 Vpk / w: twice
    // the peak, and twice what its own sinusoidal fallback (Vpk / w) returns.
    const double frequency = 100000;
    const double voltagePeakToPeak = 600;
    const double omega = 2 * std::numbers::pi * frequency;
    const double expectedPeakFluxLinkage = voltagePeakToPeak / 2 / omega;
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(
        frequency, 100e-6, 25, WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, {10});
    auto excitation = OpenMagnetics::Inputs::get_primary_excitation(inputs.get_operating_points()[0]);
    REQUIRE(excitation.get_voltage());
    REQUIRE(excitation.get_voltage()->get_waveform());

    double fromWaveform = OpenMagnetics::Inputs::calculate_max_volt_seconds(excitation);
    INFO("from waveform " << fromWaveform << " V s, Vpk / w " << expectedPeakFluxLinkage << " V s");
    CHECK_THAT(fromWaveform, Catch::Matchers::WithinRel(expectedPeakFluxLinkage, 0.01));

    SECTION("agrees with the sinusoidal fallback") {
        auto voltage = excitation.get_voltage().value();
        REQUIRE(voltage.get_processed());
        REQUIRE(voltage.get_processed()->get_peak());
        voltage.set_waveform(std::nullopt);
        auto processedOnly = excitation;
        processedOnly.set_voltage(voltage);
        double fromFallback = OpenMagnetics::Inputs::calculate_max_volt_seconds(processedOnly);
        INFO("fallback " << fromFallback << " V s");
        CHECK_THAT(fromFallback, Catch::Matchers::WithinRel(expectedPeakFluxLinkage, 0.01));
        CHECK_THAT(fromWaveform, Catch::Matchers::WithinRel(fromFallback, 0.01));
    }

    SECTION("does not depend on where the period starts, nor on a DC bias") {
        // Same 300 V sine, sampled from a third of the way into its period, riding on 5 V of DC
        // (a winding cannot hold a steady voltage; the bias is a ramp, not a flux excursion).
        const size_t numberPoints = 1024;
        std::vector<double> time(numberPoints + 1);
        std::vector<double> data(numberPoints + 1);
        for (size_t k = 0; k <= numberPoints; ++k) {
            time[k] = static_cast<double>(k) / numberPoints / frequency;
            data[k] = 5 + voltagePeakToPeak / 2 * std::sin(omega * time[k] + 2 * std::numbers::pi / 3);
        }
        Waveform waveform;
        waveform.set_time(time);
        waveform.set_data(data);
        auto voltage = excitation.get_voltage().value();
        voltage.set_waveform(waveform);
        auto shifted = excitation;
        shifted.set_voltage(voltage);
        double fromShifted = OpenMagnetics::Inputs::calculate_max_volt_seconds(shifted);
        INFO("shifted and biased " << fromShifted << " V s");
        CHECK_THAT(fromShifted, Catch::Matchers::WithinRel(expectedPeakFluxLinkage, 0.01));
    }
}
