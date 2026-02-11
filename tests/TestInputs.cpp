#include <source_location>
#include "support/CoilMesher.h"
#include "support/Painter.h"
#include "physical_models/InitialPermeability.h"
#include "processors/Inputs.h"
#include "support/Utils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
using json = nlohmann::json;
#include "TestingUtils.h"

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
    REQUIRE_THAT(excitation.get_voltage().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(692290, max_error * 692290));
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
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(746020, max_error * 746020));
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
    REQUIRE_THAT(excitation.get_current().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(465000, max_error * 465000));
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
    windingExcitation["voltage"]["waveform"]["data"] = {0, 1, 1, 0, 0};
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
        REQUIRE_THAT(excitation_primary.get_current().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(129700, max_error * 129700));
        REQUIRE_THAT(excitation_primary.get_current().value().get_processed().value().get_peak_to_peak().value(), Catch::Matchers::WithinAbs(10, max_error * 10));
        REQUIRE_THAT(excitation_primary.get_current().value().get_processed().value().get_offset(), Catch::Matchers::WithinAbs(0, max_error));

        REQUIRE_THAT(excitation_primary.get_current().value().get_harmonics().value().get_amplitudes()[0], Catch::Matchers::WithinAbs(0, max_error));
        REQUIRE_THAT(excitation_primary.get_current().value().get_harmonics().value().get_amplitudes()[1], Catch::Matchers::WithinAbs(3.833, max_error * 3.833));
        REQUIRE_THAT(excitation_primary.get_current().value().get_harmonics().value().get_frequencies()[0], Catch::Matchers::WithinAbs(0, max_error));
        REQUIRE_THAT(excitation_primary.get_current().value().get_harmonics().value().get_frequencies()[1], Catch::Matchers::WithinAbs(100000, max_error * 100000));

        REQUIRE_THAT(excitation_primary.get_voltage().value().get_processed().value().get_rms().value(), Catch::Matchers::WithinAbs(0.5, max_error * 0.5));
        REQUIRE_THAT(excitation_primary.get_voltage().value().get_processed().value().get_thd().value(), Catch::Matchers::WithinAbs(0.95, max_error * 0.92));
        REQUIRE_THAT(excitation_primary.get_voltage().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(640900, max_error * 640900));
        REQUIRE_THAT(excitation_primary.get_voltage().value().get_processed().value().get_peak_to_peak().value(), Catch::Matchers::WithinAbs(1, max_error * 1));
        REQUIRE_THAT(excitation_primary.get_voltage().value().get_processed().value().get_offset(), Catch::Matchers::WithinAbs(0, max_error * 0.25));

        REQUIRE_THAT(excitation_primary.get_voltage().value().get_harmonics().value().get_amplitudes()[0], Catch::Matchers::WithinAbs(0.2421, max_error * 0.25));
        REQUIRE_THAT(excitation_primary.get_voltage().value().get_harmonics().value().get_amplitudes()[1], Catch::Matchers::WithinAbs(0.439, max_error * 0.451));
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
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(129000, max_error * 129000));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_ac_effective_frequency().value(), Catch::Matchers::WithinAbs(129000, max_error * 129000));
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
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_effective_frequency().value(), Catch::Matchers::WithinAbs(9777, max_error * 9777));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_ac_effective_frequency().value(), Catch::Matchers::WithinAbs(129000, max_error * 129000));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_peak_to_peak().value(), Catch::Matchers::WithinAbs(1.855, max_error * 1.855));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_processed().value().get_offset(), Catch::Matchers::WithinAbs(0, 0));

    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_harmonics().value().get_amplitudes()[0], Catch::Matchers::WithinAbs(10, 10 * max_error));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_harmonics().value().get_amplitudes()[1], Catch::Matchers::WithinAbs(0.7165, max_error * 0.7165));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_harmonics().value().get_frequencies()[0], Catch::Matchers::WithinAbs(0, max_error));
    REQUIRE_THAT(excitation.get_magnetizing_current().value().get_harmonics().value().get_frequencies()[1], Catch::Matchers::WithinAbs(100000, max_error * 100000));
}

TEST_CASE("Test_Get_Duty_Cycle_Triangular", "[processor][inputs][smoke-test]") {
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

    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    auto dutyCycle = OpenMagnetics::Inputs::try_guess_duty_cycle(excitation.get_current().value().get_waveform().value());
    REQUIRE(dutyCycle == 0.25);
}

TEST_CASE("Test_Get_Duty_Cycle_Rectangular", "[processor][inputs][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["current"]["waveform"]["data"] = {7.5, 7.5, -2.5, -2.5, 7.5};
    windingExcitation["current"]["waveform"]["time"] = {0, 0.0000075, 0.0000075, 0.00001, 0.00001};
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);

    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    auto dutyCycle = OpenMagnetics::Inputs::try_guess_duty_cycle(excitation.get_current().value().get_waveform().value());
    double max_error = 0.02;

    REQUIRE_THAT(0.75, Catch::Matchers::WithinAbs(dutyCycle, max_error * 0.75));
}

TEST_CASE("Test_Get_Duty_Cycle_Custom", "[processor][inputs][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["current"]["waveform"]["data"] = {0, 3, 8, 3, 0};
    windingExcitation["current"]["waveform"]["time"] = {0, 0.0000025, 0.0000042, 0.0000075, 0.00001};
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);

    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    auto dutyCycle = OpenMagnetics::Inputs::try_guess_duty_cycle(excitation.get_current().value().get_waveform().value());

    REQUIRE(dutyCycle == 0.42);
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

TEST_CASE("Test_Try_Guess_Duty_Cycle", "[processor][inputs][smoke-test]") {
    double max_error = 0.1;
    double peakToPeak = 53.3333;
    double dutyCycle = 0.25;
    double frequency = 100000;

    Waveform voltageWaveform;
    voltageWaveform.set_time(std::vector<double>{0, dutyCycle / frequency, (dutyCycle - 0.01) / frequency, 1.0 / frequency, 1.0 / frequency});
    voltageWaveform.set_data(std::vector<double>{peakToPeak * (1 - dutyCycle), peakToPeak * (1 - dutyCycle), -peakToPeak * dutyCycle, -peakToPeak * dutyCycle, peakToPeak * (1 - dutyCycle)});

    double expectedValue = 0.25;

    auto dutyCycleGuessed = OpenMagnetics::Inputs::try_guess_duty_cycle(voltageWaveform);
    REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(dutyCycleGuessed, max_error * expectedValue));
}

TEST_CASE("Test_Try_Guess_Sinusoidal", "[processor][inputs][smoke-test]") {
    double max_error = 0.1;
    double peakToPeak = 53.3333;
    double dutyCycle = 0.5;
    double offset = 0;
    double frequency = 100000;
    auto label = WaveformLabel::SINUSOIDAL;

    Processed processed;
    processed.set_peak_to_peak(peakToPeak);
    processed.set_duty_cycle(dutyCycle);
    processed.set_offset(offset);
    processed.set_label(label);

    auto waveform = OpenMagnetics::Inputs::create_waveform(processed, frequency);
    auto guessedLabel = OpenMagnetics::Inputs::try_guess_waveform_label(waveform);

    auto calculatedProcessed = OpenMagnetics::Inputs::calculate_basic_processed_data(waveform);

    REQUIRE(magic_enum::enum_name(label) == magic_enum::enum_name(guessedLabel));
    REQUIRE_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value()));
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_duty_cycle().value(), max_error * processed.get_duty_cycle().value()));
    REQUIRE_THAT(processed.get_offset(), Catch::Matchers::WithinAbs(calculatedProcessed.get_offset(), max_error));
    REQUIRE(magic_enum::enum_name(calculatedProcessed.get_label()) == magic_enum::enum_name(processed.get_label()));
}

TEST_CASE("Test_Try_Guess_Triangular", "[processor][inputs][smoke-test]") {
    double max_error = 0.1;
    double peakToPeak = 53.3333;
    double dutyCycle = 0.25;
    double offset = 0;
    double frequency = 100000;
    auto label = WaveformLabel::TRIANGULAR;

    Processed processed;
    processed.set_peak_to_peak(peakToPeak);
    processed.set_duty_cycle(dutyCycle);
    processed.set_offset(offset);
    processed.set_label(label);

    auto waveform = OpenMagnetics::Inputs::create_waveform(processed, frequency);
    auto guessedLabel = OpenMagnetics::Inputs::try_guess_waveform_label(waveform);

    auto calculatedProcessed = OpenMagnetics::Inputs::calculate_basic_processed_data(waveform);

    REQUIRE(magic_enum::enum_name(label) == magic_enum::enum_name(guessedLabel));
    REQUIRE_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value()));
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_duty_cycle().value(), max_error * processed.get_duty_cycle().value()));
    REQUIRE_THAT(processed.get_offset(), Catch::Matchers::WithinAbs(calculatedProcessed.get_offset(), max_error));
    REQUIRE(magic_enum::enum_name(calculatedProcessed.get_label()) == magic_enum::enum_name(processed.get_label()));
}

TEST_CASE("Test_Try_Guess_Unipolar_Triangular", "[processor][inputs][smoke-test]") {
    double max_error = 0.1;
    double peakToPeak = 53.3333;
    double dutyCycle = 0.25;
    double offset = 0;
    double frequency = 100000;
    auto label = WaveformLabel::UNIPOLAR_TRIANGULAR;

    Processed processed;
    processed.set_peak_to_peak(peakToPeak);
    processed.set_duty_cycle(dutyCycle);
    processed.set_offset(offset);
    processed.set_label(label);

    auto waveform = OpenMagnetics::Inputs::create_waveform(processed, frequency);
    auto guessedLabel = OpenMagnetics::Inputs::try_guess_waveform_label(waveform);

    auto calculatedProcessed = OpenMagnetics::Inputs::calculate_processed_data(waveform, frequency);

    double expectedAverage = (peakToPeak * dutyCycle / 2);
    REQUIRE(magic_enum::enum_name(label) == magic_enum::enum_name(guessedLabel));
    REQUIRE_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value()));
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_duty_cycle().value(), max_error * processed.get_duty_cycle().value()));
    REQUIRE_THAT(0, Catch::Matchers::WithinAbs(calculatedProcessed.get_offset(), max_error));
    REQUIRE_THAT(expectedAverage, Catch::Matchers::WithinAbs(calculatedProcessed.get_average().value(), max_error * expectedAverage));
    REQUIRE(magic_enum::enum_name(calculatedProcessed.get_label()) == magic_enum::enum_name(processed.get_label()));
}

TEST_CASE("Test_Try_Guess_Unipolar_Rectangular", "[processor][inputs][smoke-test]") {
    double max_error = 0.1;
    double peakToPeak = 53.3333;
    double dutyCycle = 0.25;
    double offset = 0;
    double frequency = 100000;
    auto label = WaveformLabel::UNIPOLAR_RECTANGULAR;

    Processed processed;
    processed.set_peak_to_peak(peakToPeak);
    processed.set_duty_cycle(dutyCycle);
    processed.set_offset(offset);
    processed.set_label(label);

    auto waveform = OpenMagnetics::Inputs::create_waveform(processed, frequency);
    auto guessedLabel = OpenMagnetics::Inputs::try_guess_waveform_label(waveform);

    auto calculatedProcessed = OpenMagnetics::Inputs::calculate_processed_data(waveform, frequency);

    double expectedAverage = peakToPeak * dutyCycle;
    REQUIRE(magic_enum::enum_name(label) == magic_enum::enum_name(guessedLabel));
    REQUIRE_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value()));
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_duty_cycle().value(), max_error * processed.get_duty_cycle().value()));
    REQUIRE_THAT(0, Catch::Matchers::WithinAbs(calculatedProcessed.get_offset(), max_error));
    REQUIRE_THAT(expectedAverage, Catch::Matchers::WithinAbs(calculatedProcessed.get_average().value(), max_error * expectedAverage));
    REQUIRE(magic_enum::enum_name(calculatedProcessed.get_label()) == magic_enum::enum_name(processed.get_label()));
}

TEST_CASE("Test_Try_Guess_Rectangular", "[processor][inputs][smoke-test]") {
    double max_error = 0.1;
    double peakToPeak = 53.3333;
    double dutyCycle = 0.25;
    double offset = 0;
    double frequency = 100000;
    auto label = WaveformLabel::RECTANGULAR;

    Processed processed;
    processed.set_peak_to_peak(peakToPeak);
    processed.set_duty_cycle(dutyCycle);
    processed.set_offset(offset);
    processed.set_label(label);

    auto waveform = OpenMagnetics::Inputs::create_waveform(processed, frequency);
    auto guessedLabel = OpenMagnetics::Inputs::try_guess_waveform_label(waveform);

    auto calculatedProcessed = OpenMagnetics::Inputs::calculate_basic_processed_data(waveform);

    double expectedOffset = 0;
    REQUIRE(magic_enum::enum_name(label) == magic_enum::enum_name(guessedLabel));
    REQUIRE_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value()));
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_duty_cycle().value(), max_error * processed.get_duty_cycle().value()));
    REQUIRE_THAT(expectedOffset, Catch::Matchers::WithinAbs(calculatedProcessed.get_offset(), max_error));
    REQUIRE(magic_enum::enum_name(calculatedProcessed.get_label()) == magic_enum::enum_name(processed.get_label()));
}

TEST_CASE("Test_Try_Guess_Bipolar_Rectangular", "[processor][inputs][smoke-test]") {
    double max_error = 0.2;
    double peakToPeak = 53.3333;
    double dutyCycle = 0.5;
    double offset = 0;
    double frequency = 100000;
    auto label = WaveformLabel::BIPOLAR_RECTANGULAR;

    Processed processed;
    processed.set_peak_to_peak(peakToPeak);
    processed.set_duty_cycle(dutyCycle);
    processed.set_offset(offset);
    processed.set_label(label);

    auto waveform = OpenMagnetics::Inputs::create_waveform(processed, frequency);
    auto guessedLabel = OpenMagnetics::Inputs::try_guess_waveform_label(waveform);

    auto calculatedProcessed = OpenMagnetics::Inputs::calculate_basic_processed_data(waveform);

    double expectedOffset = 0;
    REQUIRE(magic_enum::enum_name(label) == magic_enum::enum_name(guessedLabel));
    REQUIRE_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value()));
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_duty_cycle().value(), max_error * processed.get_duty_cycle().value()));
    REQUIRE_THAT(expectedOffset, Catch::Matchers::WithinAbs(calculatedProcessed.get_offset(), max_error));
    REQUIRE(magic_enum::enum_name(calculatedProcessed.get_label()) == magic_enum::enum_name(processed.get_label()));
}

TEST_CASE("Test_Try_Guess_Bipolar_Triangular", "[processor][inputs][smoke-test]") {
    double max_error = 0.2;
    double peakToPeak = 53.3333;
    double dutyCycle = 0.5;
    double offset = 0;
    double frequency = 100000;
    auto label = WaveformLabel::BIPOLAR_TRIANGULAR;

    Processed processed;
    processed.set_peak_to_peak(peakToPeak);
    processed.set_duty_cycle(dutyCycle);
    processed.set_offset(offset);
    processed.set_label(label);

    auto waveform = OpenMagnetics::Inputs::create_waveform(processed, frequency);
    auto guessedLabel = OpenMagnetics::Inputs::try_guess_waveform_label(waveform);

    auto calculatedProcessed = OpenMagnetics::Inputs::calculate_basic_processed_data(waveform);

    double expectedOffset = 0;
    REQUIRE(magic_enum::enum_name(label) == magic_enum::enum_name(guessedLabel));
    REQUIRE_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value()));
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_duty_cycle().value(), max_error * processed.get_duty_cycle().value()));
    REQUIRE_THAT(expectedOffset, Catch::Matchers::WithinAbs(calculatedProcessed.get_offset(), max_error));
    REQUIRE(magic_enum::enum_name(calculatedProcessed.get_label()) == magic_enum::enum_name(processed.get_label()));
}

TEST_CASE("Test_Try_Guess_Flyback_Primary", "[processor][inputs][smoke-test]") {
    double max_error = 0.1;
    double peakToPeak = 50;
    double dutyCycle = 0.25;
    double offset = 10;
    double frequency = 100000;
    auto label = WaveformLabel::FLYBACK_PRIMARY;

    Processed processed;
    processed.set_peak_to_peak(peakToPeak);
    processed.set_duty_cycle(dutyCycle);
    processed.set_offset(offset);
    processed.set_label(label);

    auto waveform = OpenMagnetics::Inputs::create_waveform(processed, frequency);
    auto guessedLabel = OpenMagnetics::Inputs::try_guess_waveform_label(waveform);

    auto calculatedProcessed = OpenMagnetics::Inputs::calculate_basic_processed_data(waveform);

    REQUIRE(magic_enum::enum_name(label) == magic_enum::enum_name(guessedLabel));
    REQUIRE_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value()));
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_duty_cycle().value(), max_error * processed.get_duty_cycle().value()));
    REQUIRE(magic_enum::enum_name(calculatedProcessed.get_label()) == magic_enum::enum_name(processed.get_label()));
}

TEST_CASE("Test_Try_Guess_Flyback_Secondary", "[processor][inputs][smoke-test]") {
    double max_error = 0.1;
    double peakToPeak = 50;
    double dutyCycle = 0.25;
    double offset = 10;
    double frequency = 100000;
    auto label = WaveformLabel::FLYBACK_SECONDARY;

    Processed processed;
    processed.set_peak_to_peak(peakToPeak);
    processed.set_duty_cycle(dutyCycle);
    processed.set_offset(offset);
    processed.set_label(label);

    auto waveform = OpenMagnetics::Inputs::create_waveform(processed, frequency);
    auto guessedLabel = OpenMagnetics::Inputs::try_guess_waveform_label(waveform);

    auto calculatedProcessed = OpenMagnetics::Inputs::calculate_basic_processed_data(waveform);

    REQUIRE(magic_enum::enum_name(label) == magic_enum::enum_name(guessedLabel));
    REQUIRE_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value()));
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_duty_cycle().value(), max_error * processed.get_duty_cycle().value()));
    REQUIRE(magic_enum::enum_name(calculatedProcessed.get_label()) == magic_enum::enum_name(processed.get_label()));
}

TEST_CASE("Test_Sample_And_Compress_Sinusoidal", "[processor][inputs][smoke-test]") {
    Waveform waveform;
    std::vector<double> data;
    std::vector<double> time;
    double peakToPeak = 50;
    double offset = 0;
    size_t numberPointsSampledWaveforms = 69;
    for (size_t i = 0; i < numberPointsSampledWaveforms; ++i) {
        double angle = i * 2 * M_PI / numberPointsSampledWaveforms;
        time.push_back(angle);
        data.push_back((sin(angle) * peakToPeak / 2) + offset);
    }
    waveform.set_data(data);
    waveform.set_time(time);

    REQUIRE(!OpenMagnetics::Inputs::is_waveform_sampled(waveform));

    auto sampledWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(waveform);
    REQUIRE(OpenMagnetics::Inputs::is_waveform_sampled(sampledWaveform));
    REQUIRE(128UL == sampledWaveform.get_data().size());

    auto compressedWaveform = OpenMagnetics::Inputs::compress_waveform(sampledWaveform);
    REQUIRE(115UL == compressedWaveform.get_data().size());
    auto guessedLabel = OpenMagnetics::Inputs::try_guess_waveform_label(compressedWaveform);
    REQUIRE(magic_enum::enum_name(WaveformLabel::SINUSOIDAL) == magic_enum::enum_name(guessedLabel));
}

TEST_CASE("Test_Sample_And_Compress_Rectangular", "[processor][inputs][smoke-test]") {
    Waveform waveform;
    waveform.set_data(std::vector<double>({-2.5, 7.5, 7.5, -2.5, -2.5}));
    waveform.set_time(std::vector<double>({0, 0, 0.0000025, 0.0000025, 0.00001}));

    REQUIRE(!OpenMagnetics::Inputs::is_waveform_sampled(waveform));

    auto sampledWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(waveform);
    REQUIRE(OpenMagnetics::Inputs::is_waveform_sampled(sampledWaveform));
    REQUIRE(128UL == sampledWaveform.get_data().size());

    auto compressedWaveform = OpenMagnetics::Inputs::compress_waveform(sampledWaveform);
    REQUIRE(5UL == compressedWaveform.get_data().size());
    auto guessedLabel = OpenMagnetics::Inputs::try_guess_waveform_label(compressedWaveform);
    REQUIRE(magic_enum::enum_name(WaveformLabel::RECTANGULAR) == magic_enum::enum_name(guessedLabel));
}

TEST_CASE("Test_Sample_And_Compress_Triangular", "[processor][inputs][smoke-test]") {
    Waveform waveform;
    waveform.set_data(std::vector<double>({-5, 5, -5}));
    waveform.set_time(std::vector<double>({0, 0.0000025, 0.00001}));

    REQUIRE(!OpenMagnetics::Inputs::is_waveform_sampled(waveform));

    auto sampledWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(waveform);
    REQUIRE(OpenMagnetics::Inputs::is_waveform_sampled(sampledWaveform));
    REQUIRE(128UL == sampledWaveform.get_data().size());

    auto compressedWaveform = OpenMagnetics::Inputs::compress_waveform(sampledWaveform);
    REQUIRE(3UL == compressedWaveform.get_data().size());
    auto guessedLabel = OpenMagnetics::Inputs::try_guess_waveform_label(compressedWaveform);
    REQUIRE(magic_enum::enum_name(WaveformLabel::TRIANGULAR) == magic_enum::enum_name(guessedLabel));
}

TEST_CASE("Test_Sample_And_Compress_Unipolar_Triangular", "[processor][inputs][smoke-test]") {
    Waveform waveform;
    waveform.set_data(std::vector<double>({10, 60, 10, 10}));
    waveform.set_time(std::vector<double>({0, 0.0000025, 0.0000025, 0.00001}));

    REQUIRE(!OpenMagnetics::Inputs::is_waveform_sampled(waveform));

    auto sampledWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(waveform);
    REQUIRE(OpenMagnetics::Inputs::is_waveform_sampled(sampledWaveform));
    REQUIRE(128UL == sampledWaveform.get_data().size());

    auto compressedWaveform = OpenMagnetics::Inputs::compress_waveform(sampledWaveform);
    REQUIRE(4UL == compressedWaveform.get_data().size());
    auto guessedLabel = OpenMagnetics::Inputs::try_guess_waveform_label(compressedWaveform);
    REQUIRE(magic_enum::enum_name(WaveformLabel::UNIPOLAR_TRIANGULAR) == magic_enum::enum_name(guessedLabel));
}

TEST_CASE("Test_Sample_And_Compress_Unipolar_Rectangular", "[processor][inputs][smoke-test]") {
    Waveform waveform;
    waveform.set_data(std::vector<double>({10, 60, 60, 10, 10}));
    waveform.set_time(std::vector<double>({0, 0, 0.0000025, 0.0000025, 0.00001}));

    REQUIRE(!OpenMagnetics::Inputs::is_waveform_sampled(waveform));

    auto sampledWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(waveform);
    REQUIRE(OpenMagnetics::Inputs::is_waveform_sampled(sampledWaveform));
    REQUIRE(128UL == sampledWaveform.get_data().size());

    auto compressedWaveform = OpenMagnetics::Inputs::compress_waveform(sampledWaveform);
    REQUIRE(5UL == compressedWaveform.get_data().size());
    auto guessedLabel = OpenMagnetics::Inputs::try_guess_waveform_label(compressedWaveform);
    REQUIRE(magic_enum::enum_name(WaveformLabel::UNIPOLAR_RECTANGULAR) == magic_enum::enum_name(guessedLabel));
}

TEST_CASE("Test_Sample_And_Compress_Bipolar_Rectangular", "[processor][inputs][smoke-test]") {
    Waveform waveform;
    waveform.set_data(std::vector<double>({0, 0, 25, 25, 0, 0, -25, -25, 0, 0}));
    waveform.set_time(std::vector<double>({0, 1.25e-6, 1.25e-6, 3.75e-6, 3.75e-6, 6.25e-6, 6.25e-6, 8.75e-6, 8.75e-6, 10e-6}));

    REQUIRE(!OpenMagnetics::Inputs::is_waveform_sampled(waveform));

    auto sampledWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(waveform);
    REQUIRE(OpenMagnetics::Inputs::is_waveform_sampled(sampledWaveform));
    REQUIRE(128UL == sampledWaveform.get_data().size());

    auto compressedWaveform = OpenMagnetics::Inputs::compress_waveform(sampledWaveform);
    REQUIRE(10UL == compressedWaveform.get_data().size());
    auto guessedLabel = OpenMagnetics::Inputs::try_guess_waveform_label(compressedWaveform);
    REQUIRE(magic_enum::enum_name(WaveformLabel::BIPOLAR_RECTANGULAR) == magic_enum::enum_name(guessedLabel));
}

TEST_CASE("Test_Sample_And_Compress_Flyback_Primary", "[processor][inputs][smoke-test]") {
    Waveform waveform;
    waveform.set_data(std::vector<double>({0, 30, 80, 0, 0}));
    waveform.set_time(std::vector<double>({0, 0, 1.4e-6, 1.4e-6, 0.00001}));

    REQUIRE(!OpenMagnetics::Inputs::is_waveform_sampled(waveform));

    auto sampledWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(waveform);
    REQUIRE(OpenMagnetics::Inputs::is_waveform_sampled(sampledWaveform));
    REQUIRE(128UL == sampledWaveform.get_data().size());

    auto compressedWaveform = OpenMagnetics::Inputs::compress_waveform(sampledWaveform);
    REQUIRE(5UL == compressedWaveform.get_data().size());
    auto guessedLabel = OpenMagnetics::Inputs::try_guess_waveform_label(compressedWaveform);
    REQUIRE(magic_enum::enum_name(WaveformLabel::FLYBACK_PRIMARY) == magic_enum::enum_name(guessedLabel));
}

TEST_CASE("Test_Sample_And_Compress_Flyback_Secondary", "[processor][inputs][smoke-test]") {
    Waveform waveform;
    waveform.set_data(std::vector<double>({0, 0, 80, 30, 0}));
    waveform.set_time(std::vector<double>({0, 1.4e-6, 1.4e-6, 0.00001, 0.00001}));

    REQUIRE(!OpenMagnetics::Inputs::is_waveform_sampled(waveform));

    auto sampledWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(waveform);
    REQUIRE(OpenMagnetics::Inputs::is_waveform_sampled(sampledWaveform));
    REQUIRE(128UL == sampledWaveform.get_data().size());

    auto compressedWaveform = OpenMagnetics::Inputs::compress_waveform(sampledWaveform);
    REQUIRE(5UL == compressedWaveform.get_data().size());
    auto guessedLabel = OpenMagnetics::Inputs::try_guess_waveform_label(compressedWaveform);
    REQUIRE(magic_enum::enum_name(WaveformLabel::FLYBACK_SECONDARY) == magic_enum::enum_name(guessedLabel));
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
        double angle = i * 2 * M_PI / numberPointsSampledWaveforms;
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
    auto processed = OpenMagnetics::Inputs::calculate_processed_data(signal.get_waveform().value());
    
    auto reflectedSignal = OpenMagnetics::Inputs::reflect_waveform(signal, ratio, WaveformLabel::UNIPOLAR_RECTANGULAR);
    auto reflectedProcessed = OpenMagnetics::Inputs::calculate_processed_data(reflectedSignal.get_waveform().value());

    REQUIRE_THAT(-processed.get_average().value() * ratio, Catch::Matchers::WithinAbs(reflectedProcessed.get_average().value(), max_error * processed.get_average().value() * ratio));
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
    auto processed = OpenMagnetics::Inputs::calculate_processed_data(signal.get_waveform().value());
    
    auto reflectedSignal = OpenMagnetics::Inputs::reflect_waveform(signal, ratio, WaveformLabel::UNIPOLAR_RECTANGULAR);
    auto reflectedProcessed = OpenMagnetics::Inputs::calculate_processed_data(reflectedSignal.get_waveform().value());

    REQUIRE_THAT(-processed.get_average().value() * ratio, Catch::Matchers::WithinAbs(reflectedProcessed.get_average().value(), max_error * processed.get_average().value() * ratio));
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
    auto processed = OpenMagnetics::Inputs::calculate_processed_data(signal.get_waveform().value());
    
    auto reflectedSignal = OpenMagnetics::Inputs::reflect_waveform(signal, ratio, WaveformLabel::UNIPOLAR_TRIANGULAR);
    auto reflectedProcessed = OpenMagnetics::Inputs::calculate_processed_data(reflectedSignal.get_waveform().value());

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

    OperatingPointExcitation primaryExcitation = json::parse(R"({"current":{"harmonics":{"amplitudes":[4.558291965060757,6.529972148854835,1.789145464052662,1.942710261723393,1.1154397290543905,1.0065940903149284,0.8988582408501911,0.5788895059301848,0.7630891566656498,0.35353473619244935,0.6456068660055856,0.2676407747861538,0.5321600974207273,0.27460271311591355,0.4212435511767113,0.3071864237685926,0.31636527222157335,0.33025480179543903,0.22484933866699708,0.3342584403228846,0.1602304192861384,0.31871495949628625,0.1407347185986853,0.2864343451372125,0.16092094387433578,0.24187051785845007,0.19181915742756442,0.19085711275346715,0.21571905868290692,0.141477779475246,0.2263289212948745,0.10671871649160397,0.22229630688335306,0.10207527130577318,0.20470089894984397,0.12300075015911581,0.17626009130576112,0.14989831972782144,0.1412892563096207,0.17104729068836264,0.10662686582719509,0.1817193278103439,0.084043934810831,0.1804622345981364,0.08623921969773558,0.16772355314774284,0.10740317057400037,0.14545650387081813,0.13211333071996756,0.11727666410506915,0.15153636839094595,0.0895735293816873,0.16177942374108412,0.07382571026303117,0.16140346655602664,0.08062965411618814,0.1505782829167501,0.10267745807614781,0.1309155701859174,0.12680055354222536,0.1058030340183861,0.14569473884415238,0.08179428173758013,0.15587714970397382],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]},"processed":{"acEffectiveFrequency":621718.3722552633,"average":5.42279561360676,"dutyCycle":0.4609374999999997,"effectiveFrequency":531425.8014807898,"label":"Flyback Primary","offset":5.448133680555522,"peak":14.514900394303213,"peakToPeak":9.223090277777805,"rms":6.997011716385439,"thd":0.5671816768565314},"waveform":{"ancillaryLabel":null,"data":[0,5.448133680555522,14.671223958333327,0,0],"numberPeriods":null,"time":[0,0,0.000004609374999999998,0.000004609374999999998,0.00001]}},"frequency":100000,"magneticFieldStrength":null,"magneticFluxDensity":null,"magnetizingCurrent":null,"name":"Primary winding excitation","voltage":{"harmonics":{"amplitudes":[9.22309027777781,12.643860752984521,1.5482098237509736,3.96519064048262,1.5036230954105976,2.088180163711471,1.4309295961552153,1.194492029523459,1.332493501826718,0.6415741456288709,1.2115023045865507,0.2568357589566538,1.0718500529421924,0.024460284513222373,0.9179949822505038,0.23169804853214843,0.7547969368000894,0.3802751568040793,0.5873407469878427,0.47948718542917224,0.4207522514685872,0.5358937098386783,0.2600139243728212,0.5548945409090823,0.10978707149082016,0.5414997594159607,0.025752703324263687,0.5006963718234527,0.1430605791849741,0.4375948229458343,0.2393427194140055,0.35744706876955396,0.312647128060265,0.2655866960528321,0.3619218291128731,0.16732178127549502,0.38703829474775314,0.06780108548165185,0.3887793859782987,0.028131331656939335,0.36879243284993984,0.11608129419643629,0.3295094194096947,0.19223993817404458,0.27403749968881713,0.2535015787050763,0.2060242064300918,0.2975577510911559,0.12950268067955925,0.3229626417701206,0.048723010567191104,0.32916660187526037,0.03202370751910006,0.3165159206843967,0.10858170269138588,0.2862185126841811,0.17710037177346274,0.24027660937080772,0.23419100538841045,0.18138890635478225,0.2770642907631915,0.11282585612072807,0.3036432370169665,0.03828287613148918],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]},"processed":{"acEffectiveFrequency":593227.8773244782,"average":9.223090277777818,"dutyCycle":0.46,"effectiveFrequency":496491.712181633,"label":"Custom","offset":10.004708097928473,"peak":20.009416195856947,"peakToPeak":20.009416195856947,"rms":13.584868493291284,"thd":0.49423464655132765},"waveform":{"ancillaryLabel":null,"data":[20.009416195856947,20.009416195856947,0,0],"numberPeriods":null,"time":[0,0.000004609374999999997,0.000004609374999999997,0.000009999999999999999]}}})");

    OperatingPointExcitation excitationOfThisWinding(primaryExcitation);
    auto currentSignalDescriptorProcessed = OpenMagnetics::Inputs::calculate_basic_processed_data(primaryExcitation.get_current().value().get_waveform().value());
    auto reflectedCurrentSignalDescriptor = OpenMagnetics::Inputs::reflect_waveform(primaryExcitation.get_current().value(), turnRatio, currentSignalDescriptorProcessed.get_label());
    auto processed = OpenMagnetics::Inputs::calculate_basic_processed_data(reflectedCurrentSignalDescriptor.get_waveform().value());
    REQUIRE(processed.get_label() == WaveformLabel::FLYBACK_SECONDARY);
}

TEST_CASE("Test_Flyback_Json", "[processor][inputs][smoke-test]") {
    std::string masString = R"({"inputs": {"designRequirements": {"insulation": {"altitude": {"maximum": 1000 }, "cti": "Group IIIB", "pollutionDegree": "P2", "overvoltageCategory": "OVC-II", "insulationType": "Basic", "mainSupplyVoltage": {"nominal": 120, "minimum": 20, "maximum": 450 }, "standards": ["IEC 60664-1"] }, "leakageInductance": [{"maximum": 0.00000135 }, {"maximum": 0.00000135 } ], "magnetizingInductance": {"maximum": 0.000088, "minimum": 0.000036, "nominal": 0.000039999999999999996 }, "market": "Commercial", "maximumDimensions": {"width": null, "height": 0.02, "depth": null }, "maximumWeight": null, "name": "My Design Requirements", "operatingTemperature": null, "strayCapacitance": [{"maximum": 5e-11 }, {"maximum": 5e-11 } ], "terminalType": ["Pin", "Pin", "Pin"], "topology": "Flyback Converter", "turnsRatios": [{"nominal": 1.2 }, {"nominal": 1.2 } ] }, "operatingPoints": [{"name": "Operating Point No. 1", "conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"name": "Primary winding excitation", "frequency": 36000, "current": {"waveform": {"ancillaryLabel": null, "data": [0, 0, 5.5, 0, 0 ], "numberPeriods": null, "time": [0, 0, 0.000011111111111111112, 0.000011111111111111112, 0.00002777777777777778 ] }, "processed": {"dutyCycle": 0.4, "peakToPeak": 5.5, "offset": 0, "label": "Flyback Primary", "acEffectiveFrequency": 247671.20334838802, "effectiveFrequency": 224591.09893738205, "peak": 5.478515625000006, "rms": 2.025897603604934, "thd": 0.8123116860415294 }, "harmonics": {"amplitudes": [1.1128234863281259, 1.8581203560604207, 1.051281221517565, 0.5159082238421047, 0.4749772671818343, 0.35101798151542746, 0.28680520988943464, 0.26606329358167974, 0.2109058168232097, 0.2038324242201338, 0.17683391011628416, 0.15977271740890614, 0.15330090420784062, 0.13359870506402888, 0.1311678604251054, 0.11938227872282572, 0.11196939263234992, 0.10881687126769939, 0.09881440376099933, 0.09781517962337143, 0.09113803697574208, 0.08720137281111257, 0.08534603089826898, 0.07933100843010468, 0.07894777491614259, 0.0746052794656816, 0.07229734773884018, 0.07109041486606842, 0.0670974294019412, 0.06703716258861198, 0.06396085871081211, 0.0625474905371643, 0.06171647321147238, 0.05889230550248535, 0.059028047581821314, 0.05671822612277054, 0.05585484049699437, 0.055265364667603324, 0.0531817689812044, 0.053452567468547284, 0.05164394009764037, 0.051150939884780554, 0.05072953903613357, 0.049149522805038266, 0.049525208170729855, 0.04806384343218855, 0.04784167797419847, 0.047545488047249944, 0.04632889186784806, 0.046795775944523485, 0.045586400651502844, 0.045578916197115495, 0.045382141359988915, 0.04444264358178238, 0.04499725549277023, 0.04397809673849301, 0.04415534166549457, 0.044042673430388084, 0.04332652563384539, 0.04397294643650677, 0.043102730085842, 0.043452676225114936, 0.04341597260904029, 0.04289100286968138 ], "frequencies": [0, 36000, 72000, 108000, 144000, 180000, 216000, 252000, 288000, 324000, 360000, 396000, 432000, 468000, 504000, 540000, 576000, 612000, 648000, 684000, 720000, 756000, 792000, 828000, 864000, 900000, 936000, 972000, 1008000, 1044000, 1080000, 1116000, 1152000, 1188000, 1224000, 1260000, 1296000, 1332000, 1368000, 1404000, 1440000, 1476000, 1512000, 1548000, 1584000, 1620000, 1656000, 1692000, 1728000, 1764000, 1800000, 1836000, 1872000, 1908000, 1944000, 1980000, 2016000, 2052000, 2088000, 2124000, 2160000, 2196000, 2232000, 2268000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-8, 12, 12, -8, -8 ], "numberPeriods": null, "time": [0, 0, 0.000011111111111111112, 0.000011111111111111112, 0.00002777777777777778 ] }, "processed": {"dutyCycle": 0.4, "peakToPeak": 20, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 217445.61780293446, "effectiveFrequency": 217445.06394458748, "peak": 12, "rms": 9.79157801378307, "thd": 0.5579294461931013 }, "harmonics": {"amplitudes": [0.03125, 12.090982167347164, 3.7938629699811433, 2.4460154592542187, 3.050934294867692, 0.06265085868474476, 2.0052584267980813, 1.1245773748206938, 0.8899247078195547, 1.3746157728063957, 0.06310658027096225, 1.0931619934351904, 0.682943044655975, 0.5329821111108219, 0.8998037891938968, 0.06387675236155328, 0.7544417382415929, 0.5036912732460003, 0.37575767434186663, 0.6781566555968945, 0.06497787099509279, 0.579613707579126, 0.4082105683349906, 0.28790910341678466, 0.5516771452281819, 0.06643416431010025, 0.4742269329016237, 0.35024440150573616, 0.23220844850577863, 0.4713901788294872, 0.0682788501833682, 0.40480954492681936, 0.3125, 0.19399945657286152, 0.41719035005974064, 0.07055595095604376, 0.35652867811806305, 0.28709106236926035, 0.16634675530328177, 0.3793488450042862, 0.07332285477181877, 0.32183350571658953, 0.2699539504381825, 0.14554626367323475, 0.3526337777044684, 0.0766539126369276, 0.29650814745710025, 0.25883312123614344, 0.12944173824159266, 0.33403609127672834, 0.08064551519955634, 0.27805045130481554, 0.25243555700995646, 0.1166943265232067, 0.3217665194149353, 0.08542334196565987, 0.2649247078195551, 0.250039428516803, 0.10643002348369442, 0.31477312120987916, 0.09115288771096393, 0.25618941041585475, 0.2513050615401854, 0.098055076230366 ], "frequencies": [0, 36000, 72000, 108000, 144000, 180000, 216000, 252000, 288000, 324000, 360000, 396000, 432000, 468000, 504000, 540000, 576000, 612000, 648000, 684000, 720000, 756000, 792000, 828000, 864000, 900000, 936000, 972000, 1008000, 1044000, 1080000, 1116000, 1152000, 1188000, 1224000, 1260000, 1296000, 1332000, 1368000, 1404000, 1440000, 1476000, 1512000, 1548000, 1584000, 1620000, 1656000, 1692000, 1728000, 1764000, 1800000, 1836000, 1872000, 1908000, 1944000, 1980000, 2016000, 2052000, 2088000, 2124000, 2160000, 2196000, 2232000, 2268000 ] } } }, {"name": "Primary winding excitation", "frequency": 36000, "current": {"waveform": {"ancillaryLabel": null, "data": [0, 0, 6.6, 0, 0 ], "numberPeriods": null, "time": [0, 0.00001388888888888889, 0.00001388888888888889, 0.00002777777777777778, 0.00002777777777777778 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 6.6, "offset": 0, "label": "Flyback Secondary", "acEffectiveFrequency": 234518.76473649617, "effectiveFrequency": 205388.02000937346, "peak": 6.599999999999992, "rms": 2.7259938518765505, "thd": 0.6765065136547781 }, "harmonics": {"amplitudes": [1.6757812500000036, 2.5183033006556617, 1.0508445877408874, 0.7271846526610849, 0.5260559513023173, 0.43112241679715924, 0.35140961815458516, 0.30751649536118214, 0.2643006555483416, 0.23957382058929272, 0.21220866665144478, 0.19664985565967041, 0.17762735700582763, 0.16712773941846806, 0.15305447605739808, 0.14562523024961643, 0.13473930575287565, 0.12930491961017093, 0.12059849313526576, 0.11652868809047452, 0.10938238636775419, 0.10628456966336206, 0.10029609952951436, 0.0979137380509861, 0.09281004801094254, 0.09096893358288156, 0.08655793218347599, 0.08513597099573951, 0.08127844554340995, 0.08018779841976777, 0.07678030294711186, 0.07595637558884157, 0.07292038680986246, 0.07231479790644763, 0.06958961021254718, 0.06916556182624073, 0.06670349640899874, 0.06643264501653386, 0.06419573768805274, 0.06405603077180153, 0.062013691465176275, 0.0619878418001075, 0.060115168596257366, 0.060189559664044844, 0.058466103590436726, 0.05862999271489314, 0.05703883927855185, 0.05728377042425536, 0.05581084782757627, 0.05613021480841572, 0.054763767206133816, 0.05515248676284005, 0.05388266966565528, 0.05433693628512042, 0.05315550383577058, 0.05367260658830061, 0.05257266909511635, 0.05315085657150875, 0.052126692739434176, 0.05276507626997662, 0.05181198888815346, 0.05251047719660189, 0.051624684193004804, 0.05238394486208969 ], "frequencies": [0, 36000, 72000, 108000, 144000, 180000, 216000, 252000, 288000, 324000, 360000, 396000, 432000, 468000, 504000, 540000, 576000, 612000, 648000, 684000, 720000, 756000, 792000, 828000, 864000, 900000, 936000, 972000, 1008000, 1044000, 1080000, 1116000, 1152000, 1188000, 1224000, 1260000, 1296000, 1332000, 1368000, 1404000, 1440000, 1476000, 1512000, 1548000, 1584000, 1620000, 1656000, 1692000, 1728000, 1764000, 1800000, 1836000, 1872000, 1908000, 1944000, 1980000, 2016000, 2052000, 2088000, 2124000, 2160000, 2196000, 2232000, 2268000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-8.335, 8.335, 8.335, -8.335, -8.335 ], "numberPeriods": null, "time": [0, 0, 0.00001388888888888889, 0.00001388888888888889, 0.00002777777777777778 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 16.67, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 212934.79296426187, "effectiveFrequency": 212921.79129775826, "peak": 8.335, "rms": 8.33500000000001, "thd": 0.48331514845248513 }, "harmonics": {"amplitudes": [0.130234375, 10.610320564806702, 0.26046875, 3.531088691718988, 0.26046875, 2.1118248335514735, 0.26046875, 1.5011183941457937, 0.26046875, 1.1599200168904427, 0.26046875, 0.9412131219905466, 0.26046875, 0.7884501596246478, 0.26046875, 0.6752398432920929, 0.26046875, 0.587608084557174, 0.26046875, 0.5174625323181874, 0.26046875, 0.4597916456948479, 0.26046875, 0.4113251048108877, 0.26046875, 0.36983720876947807, 0.26046875, 0.333759858425148, 0.26046875, 0.3019556674257579, 0.26046875, 0.2735788835440449, 0.26046875, 0.24798686524225744, 0.26046875, 0.2246818889174983, 0.26046875, 0.20327180760048874, 0.26046875, 0.18344279082219334, 0.26046875, 0.16494001687000165, 0.26046875, 0.14755372430492025, 0.26046875, 0.1311089508695984, 0.26046875, 0.1154578561962637, 0.26046875, 0.10047388405843227, 0.26046875, 0.0860472521926578, 0.26046875, 0.0720814108318868, 0.26046875, 0.058490213754947895, 0.26046875, 0.04519561547653339, 0.26046875, 0.032125756193742694, 0.26046875, 0.019213329273115587, 0.26046875, 0.006394148914942832 ], "frequencies": [0, 36000, 72000, 108000, 144000, 180000, 216000, 252000, 288000, 324000, 360000, 396000, 432000, 468000, 504000, 540000, 576000, 612000, 648000, 684000, 720000, 756000, 792000, 828000, 864000, 900000, 936000, 972000, 1008000, 1044000, 1080000, 1116000, 1152000, 1188000, 1224000, 1260000, 1296000, 1332000, 1368000, 1404000, 1440000, 1476000, 1512000, 1548000, 1584000, 1620000, 1656000, 1692000, 1728000, 1764000, 1800000, 1836000, 1872000, 1908000, 1944000, 1980000, 2016000, 2052000, 2088000, 2124000, 2160000, 2196000, 2232000, 2268000 ] } } }, {"name": "Primary winding excitation", "frequency": 36000, "current": {"waveform": {"ancillaryLabel": null, "data": [0, 0, 0.08, 0, 0 ], "numberPeriods": null, "time": [0, 0.00001388888888888889, 0.00001388888888888889, 0.00002777777777777778, 0.00002777777777777778 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 0.08, "offset": 0, "label": "Flyback Secondary", "acEffectiveFrequency": 234518.7647364962, "effectiveFrequency": 205388.02000937346, "peak": 0.0799999999999999, "rms": 0.033042349719715765, "thd": 0.6765065136547783 }, "harmonics": {"amplitudes": [0.020312500000000046, 0.030524888492795898, 0.012737510154435002, 0.008814359426194969, 0.006376435773361423, 0.005225726264207993, 0.004259510523085882, 0.003727472671044633, 0.0032036443096768678, 0.0029039250980520334, 0.002572226262441755, 0.002383634614056611, 0.0021530588727979106, 0.0020257907808299163, 0.001855205770392704, 0.0017651543060559591, 0.0016332037060954637, 0.0015673323589111632, 0.0014617999167910998, 0.0014124689465512055, 0.0013258471074879294, 0.0012882978141013591, 0.0012157102973274472, 0.0011868331884968021, 0.0011249702789205155, 0.0011026537403985642, 0.0010491870567694056, 0.0010319511635847216, 0.0009851932793140603, 0.0009719733141790018, 0.0009306703387528721, 0.0009206833404708108, 0.0008838834764831795, 0.0008765430049266376, 0.0008435104268187531, 0.0008383704463786751, 0.0008085272291999843, 0.0008052441820185928, 0.0007781301537945793, 0.0007764367366278984, 0.0007516811086688035, 0.0007513677793952425, 0.0007286687102576648, 0.0007295704201702394, 0.0007086800435204447, 0.0007106665783623409, 0.0006913798700430537, 0.0006943487324152148, 0.0006764951251827443, 0.0006803662401020098, 0.0006638032388622283, 0.0006685149910647278, 0.0006531232686746098, 0.0006586295307287323, 0.0006443091374032788, 0.0006505770495551569, 0.0006372444738801973, 0.0006442528069273779, 0.0006318386998719287, 0.0006395766820603229, 0.0006280241077351938, 0.0006364906326860843, 0.0006257537477939962, 0.0006349569074192687 ], "frequencies": [0, 36000, 72000, 108000, 144000, 180000, 216000, 252000, 288000, 324000, 360000, 396000, 432000, 468000, 504000, 540000, 576000, 612000, 648000, 684000, 720000, 756000, 792000, 828000, 864000, 900000, 936000, 972000, 1008000, 1044000, 1080000, 1116000, 1152000, 1188000, 1224000, 1260000, 1296000, 1332000, 1368000, 1404000, 1440000, 1476000, 1512000, 1548000, 1584000, 1620000, 1656000, 1692000, 1728000, 1764000, 1800000, 1836000, 1872000, 1908000, 1944000, 1980000, 2016000, 2052000, 2088000, 2124000, 2160000, 2196000, 2232000, 2268000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-8.335, 8.335, 8.335, -8.335, -8.335 ], "numberPeriods": null, "time": [0, 0, 0.00001388888888888889, 0.00001388888888888889, 0.00002777777777777778 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 16.67, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 212934.79296426187, "effectiveFrequency": 212921.79129775826, "peak": 8.335, "rms": 8.33500000000001, "thd": 0.48331514845248513 }, "harmonics": {"amplitudes": [0.130234375, 10.610320564806702, 0.26046875, 3.531088691718988, 0.26046875, 2.1118248335514735, 0.26046875, 1.5011183941457937, 0.26046875, 1.1599200168904427, 0.26046875, 0.9412131219905466, 0.26046875, 0.7884501596246478, 0.26046875, 0.6752398432920929, 0.26046875, 0.587608084557174, 0.26046875, 0.5174625323181874, 0.26046875, 0.4597916456948479, 0.26046875, 0.4113251048108877, 0.26046875, 0.36983720876947807, 0.26046875, 0.333759858425148, 0.26046875, 0.3019556674257579, 0.26046875, 0.2735788835440449, 0.26046875, 0.24798686524225744, 0.26046875, 0.2246818889174983, 0.26046875, 0.20327180760048874, 0.26046875, 0.18344279082219334, 0.26046875, 0.16494001687000165, 0.26046875, 0.14755372430492025, 0.26046875, 0.1311089508695984, 0.26046875, 0.1154578561962637, 0.26046875, 0.10047388405843227, 0.26046875, 0.0860472521926578, 0.26046875, 0.0720814108318868, 0.26046875, 0.058490213754947895, 0.26046875, 0.04519561547653339, 0.26046875, 0.032125756193742694, 0.26046875, 0.019213329273115587, 0.26046875, 0.006394148914942832 ], "frequencies": [0, 36000, 72000, 108000, 144000, 180000, 216000, 252000, 288000, 324000, 360000, 396000, 432000, 468000, 504000, 540000, 576000, 612000, 648000, 684000, 720000, 756000, 792000, 828000, 864000, 900000, 936000, 972000, 1008000, 1044000, 1080000, 1116000, 1152000, 1188000, 1224000, 1260000, 1296000, 1332000, 1368000, 1404000, 1440000, 1476000, 1512000, 1548000, 1584000, 1620000, 1656000, 1692000, 1728000, 1764000, 1800000, 1836000, 1872000, 1908000, 1944000, 1980000, 2016000, 2052000, 2088000, 2124000, 2160000, 2196000, 2232000, 2268000 ] } } } ] } ] }, "magnetic": {"coil": {"functionalDescription": [{"name": "Primary"}, {"name": "Secondary"}, {"name": "Tertiary"} ] } } })";
    json masJson = json::parse(masString);

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
        REQUIRE(16384U == excitation.get_current()->get_harmonics().value().get_frequencies().size());
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

    #ifdef ENABLE_MATPLOTPP
        auto outFile = outputFilePath;
        outFile.append("Test_Standardize_Waveform.svg");
        std::filesystem::remove(outFile);   
        Painter painter(outFile, false, true);
        painter.paint_waveform(standardSignalDescriptor.get_waveform().value());
        painter.export_svg();
    #endif

}
}  // namespace
