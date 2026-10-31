#include "RandomUtils.h"
#include <source_location>
#include "support/Painter.h"
#include "support/Utils.h"
#include "processors/MagneticSimulator.h"
#include "physical_models/Temperature.h"
#include "physical_models/InitialPermeability.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/Reluctance.h"
#include "advisers/MagneticFilter.h"
#include "advisers/MagneticAdviser.h"
#include "processors/Inputs.h"
#include "processors/ResultsAnalyzer.h"
#include "support/Exceptions.h"
#include "TestingUtils.h"
#include "Fixtures.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_predicate.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <cmath>
#include <set>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

namespace {
    auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");

    // Shared sanity assertions for a simulated MAS: one output per operating point, and every
    // output must carry a finite positive magnetizing inductance, core losses and winding
    // losses. Broad physical bounds only — no goldens (a bench magnetic dissipating 10 kW is
    // as wrong as a NaN).
    void check_simulation_outputs_sane(OpenMagnetics::Mas& mas) {
        REQUIRE(mas.get_outputs().size() > 0);
        REQUIRE(mas.get_outputs().size() == mas.get_inputs().get_operating_points().size());
        for (size_t opIdx = 0; opIdx < mas.get_outputs().size(); ++opIdx) {
            INFO("Operating point " << opIdx);
            auto& output = mas.get_outputs()[opIdx];

            REQUIRE(output.get_inductance());
            double magnetizingInductance = resolve_dimensional_values(
                output.get_inductance()->get_magnetizing_inductance().get_magnetizing_inductance());
            CHECK(std::isfinite(magnetizingInductance));
            CHECK(magnetizingInductance > 0);

            REQUIRE(output.get_core_losses());
            double coreLosses = output.get_core_losses()->get_core_losses();
            CHECK(std::isfinite(coreLosses));
            CHECK(coreLosses > 0);
            CHECK(coreLosses < 1e4);

            REQUIRE(output.get_winding_losses());
            double windingLosses = output.get_winding_losses()->get_winding_losses();
            CHECK(std::isfinite(windingLosses));
            CHECK(windingLosses > 0);
            CHECK(windingLosses < 1e4);
        }
    }

    TEST_CASE("MagneticSimulator", "[processor][magnetic-simulator]") {

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {24, 24, 48};
        std::vector<int64_t> numberParallels = {1, 1, 1};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 57026;
        double magnetizingInductance = 10e-6;
        double temperature = 25;
        WaveformLabel waveShape = WaveformLabel::TRIANGULAR;
        double peakToPeak = 10;
        double dutyCycle = 0.5;
        double dcCurrent = 0;


        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        OpenMagnetics::Mas masMagnetic;
        inputs.process();

        MagneticAdviser MagneticAdviser;
        // MagneticAdviser.set_interleaving_level(interleavingLevel);
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 1);
        REQUIRE(masMagnetics.size() > 0);
    }

    TEST_CASE("MagneticSimulatorJsonHV", "[processor][magnetic-simulator]") {
        std::string masString = R"({"inputs": {"designRequirements": {"insulation": {"altitude": {"maximum": 1000 }, "cti": "groupIIIB", "pollutionDegree": "PD2", "overvoltageCategory": "II", "insulationType": "basic", "mainSupplyVoltage": {"nominal": null, "minimum": 100, "maximum": 815 }, "standards": ["IEC 60664-1"] }, "leakageInductance": [{"maximum": 0.0000027 }, {"maximum": 7.2e-7 } ], "magnetizingInductance": {"maximum": null, "minimum": null, "nominal": 0.0007 }, "market": "commercial", "maximumDimensions": {"width": null, "height": 0.02, "depth": null }, "maximumWeight": null, "name": "My Design Requirements", "operatingTemperature": null, "strayCapacitance": [{"maximum": 5e-11 }, {"maximum": 5e-11 } ], "terminalType": ["pin", "pin", "pin"], "topology": "flybackConverter", "turnsRatios": [{"nominal": 3.6 }, {"nominal": 7 } ] }, "operatingPoints": [{"name": "Operating Point No. 1", "conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"name": "Primary winding excitation", "frequency": 42000, "current": {"waveform": {"ancillaryLabel": null, "data": [0, 0, 1.1, 0, 0 ], "numberPeriods": null, "time": [0, 0, 0.000007, 0.000007, 0.00002380952380952381 ] }, "processed": {"dutyCycle": 0.294, "peakToPeak": 1.1, "offset": 0, "label": "flybackPrimary", "acEffectiveFrequency": 320324.70197566116, "effectiveFrequency": 299906.6380380921, "peak": 1.0815263605442182, "rms": 0.34251379452955516, "thd": 1.070755363168311 }, "harmonics": {"amplitudes": [0.16053906914328236, 0.29202196202762337, 0.21752279286410783, 0.1310845319633303, 0.07713860380309742, 0.06991594662329052, 0.06385479128189217, 0.04886550825820873, 0.041656942163512685, 0.04106933856014564, 0.036071164493133245, 0.03055299653431604, 0.02975135660550652, 0.028391215457215468, 0.02482068925428346, 0.023301839916889198, 0.023120564125526488, 0.021217296242203084, 0.019422484452150467, 0.019302597014669035, 0.018564086064230358, 0.016990298638197066, 0.016530128240448857, 0.016408746965119517, 0.015342013775419296, 0.014571250724490458, 0.014600261585603161, 0.014080719038131833, 0.013228104677648548, 0.013112646237623088, 0.012998033401152062, 0.012300975273665412, 0.011953654437324483, 0.012016309039770377, 0.011608286047010533, 0.011111174591424324, 0.011137874398755002, 0.01101797292674616, 0.010530630511366074, 0.01040049906973778, 0.01046208177619542, 0.010124703775041825, 0.009838808475080322, 0.009931053977331523, 0.0098034018119726, 0.009456464631666566, 0.009454698580407995, 0.009502890216301076, 0.009218272478811467, 0.00907613407579727, 0.009199516076645192, 0.00906471876755199, 0.008824889331959554, 0.008907818118804589, 0.008937898630441313, 0.00869876368067147, 0.008666288051213203, 0.008803433740422266, 0.008663604885115571, 0.008515486625290338, 0.008660125560530002, 0.008670306514975606, 0.008475927879718145, 0.008536460912837853 ], "frequencies": [0, 42000, 84000, 126000, 168000, 210000, 252000, 294000, 336000, 378000, 420000, 462000, 504000, 546000, 588000, 630000, 672000, 714000, 756000, 798000, 840000, 882000, 924000, 966000, 1008000, 1050000, 1092000, 1134000, 1176000, 1218000, 1260000, 1302000, 1344000, 1386000, 1428000, 1470000, 1512000, 1554000, 1596000, 1638000, 1680000, 1722000, 1764000, 1806000, 1848000, 1890000, 1932000, 1974000, 2016000, 2058000, 2100000, 2142000, 2184000, 2226000, 2268000, 2310000, 2352000, 2394000, 2436000, 2478000, 2520000, 2562000, 2604000, 2646000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-29.4, 70.6, 70.6, -29.4, -29.4 ], "numberPeriods": null, "time": [0, 0, 0.000007, 0.000007, 0.00002380952380952381 ] }, "processed": {"dutyCycle": 0.294, "peakToPeak": 100, "offset": 0, "label": "rectangular", "acEffectiveFrequency": 273993.1760953399, "effectiveFrequency": 273985.04815560044, "peak": 70.6, "rms": 45.33538904652732, "thd": 0.7943280084778311 }, "harmonics": {"amplitudes": [0.4937500000000009, 50.19273146886441, 30.88945734421223, 8.60726309175931, 7.514577540357868, 12.576487416046705, 7.890222788854489, 0.6723406687769768, 6.659332230305735, 6.771471403422372, 2.1663923211140617, 3.134209960512532, 5.356728264586365, 3.4347141143172375, 0.6805380572551197, 3.777520625726094, 3.7722086912079633, 1.0283716619016732, 2.1769829843489013, 3.474177141872247, 2.1027712520453545, 0.6945565902694819, 2.7474727308098412, 2.608727438661375, 0.5486770349329464, 1.7724349572072589, 2.6198081625264167, 1.4623766283543616, 0.71496628279737, 2.2318766983434566, 1.995656440067472, 0.2773827898281461, 1.562500000000001, 2.1411881515502134, 1.0841275629577949, 0.7426431917434089, 1.9342807953981735, 1.620450382955772, 0.09545259242683589, 1.4471144811333827, 1.8430944173488684, 0.8315259963159228, 0.7788650462348548, 1.7521633971823452, 1.3695427043271637, 0.042929496561759636, 1.388303796831737, 1.6472150875774871, 0.6472086912079588, 0.8254655623504514, 1.6415464278026473, 1.1917877639027972, 0.16004318769545572, 1.3696996040394052, 1.5166136111631239, 0.5023333497836172, 0.8850851521101957, 1.581550490460967, 1.0607923422947243, 0.2691521002455623, 1.3846695378376832, 1.4323337135341552, 0.38011439553529425, 0.9615889729028446 ], "frequencies": [0, 42000, 84000, 126000, 168000, 210000, 252000, 294000, 336000, 378000, 420000, 462000, 504000, 546000, 588000, 630000, 672000, 714000, 756000, 798000, 840000, 882000, 924000, 966000, 1008000, 1050000, 1092000, 1134000, 1176000, 1218000, 1260000, 1302000, 1344000, 1386000, 1428000, 1470000, 1512000, 1554000, 1596000, 1638000, 1680000, 1722000, 1764000, 1806000, 1848000, 1890000, 1932000, 1974000, 2016000, 2058000, 2100000, 2142000, 2184000, 2226000, 2268000, 2310000, 2352000, 2394000, 2436000, 2478000, 2520000, 2562000, 2604000, 2646000 ] } } }, {"name": "Primary winding excitation", "frequency": 42000, "current": {"waveform": {"ancillaryLabel": null, "data": [0, 0, 3.96, 0, 0 ], "numberPeriods": null, "time": [0, 0.000011904761904761905, 0.000011904761904761905, 0.00002380952380952381, 0.00002380952380952381 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 3.96, "offset": 0, "label": "flybackSecondary", "acEffectiveFrequency": 273605.22552591236, "effectiveFrequency": 239619.35667760254, "peak": 3.9599999999999986, "rms": 1.635596311125932, "thd": 0.6765065136547783 }, "harmonics": {"amplitudes": [1.0054687500000032, 1.5109819803933975, 0.6305067526445318, 0.43631079159665215, 0.3156335707813906, 0.25867345007829645, 0.2108457708927511, 0.1845098972167098, 0.15858039332900484, 0.14374429235357594, 0.12732519999086678, 0.11798991339580255, 0.1065764142034965, 0.10027664365108102, 0.09183268563443873, 0.08737513814977015, 0.08084358345172529, 0.0775829517661027, 0.07235909588115942, 0.06991721285428486, 0.06562943182065245, 0.06377074179801741, 0.060177659717708575, 0.05874824283059186, 0.05568602880656547, 0.054581360149729106, 0.051934759310085545, 0.0510815825974439, 0.0487670673260459, 0.048112679051860716, 0.04606818176826717, 0.04557382535330518, 0.043752232085917345, 0.04338887874386873, 0.04175376612752825, 0.04149933709574444, 0.04002209784539918, 0.03985958700992044, 0.038517442612831645, 0.03843361846308103, 0.03720821487910577, 0.037192705080064635, 0.03606910115775439, 0.036113735798427, 0.03507966215426201, 0.03517799562893598, 0.03422330356713118, 0.0343702622545533, 0.03348650869654576, 0.03367812888504953, 0.03285826032368028, 0.03309149205770409, 0.032329601799393176, 0.032602161771072315, 0.0318933023014623, 0.032203563952980456, 0.031543601457069734, 0.03189051394290526, 0.03127601564366044, 0.03165904576198601, 0.03108719333289207, 0.0315062863179612, 0.03097481051580278, 0.03143036691725387 ], "frequencies": [0, 42000, 84000, 126000, 168000, 210000, 252000, 294000, 336000, 378000, 420000, 462000, 504000, 546000, 588000, 630000, 672000, 714000, 756000, 798000, 840000, 882000, 924000, 966000, 1008000, 1050000, 1092000, 1134000, 1176000, 1218000, 1260000, 1302000, 1344000, 1386000, 1428000, 1470000, 1512000, 1554000, 1596000, 1638000, 1680000, 1722000, 1764000, 1806000, 1848000, 1890000, 1932000, 1974000, 2016000, 2058000, 2100000, 2142000, 2184000, 2226000, 2268000, 2310000, 2352000, 2394000, 2436000, 2478000, 2520000, 2562000, 2604000, 2646000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-13.89, 13.89, 13.89, -13.89, -13.89 ], "numberPeriods": null, "time": [0, 0, 0.000011904761904761905, 0.000011904761904761905, 0.00002380952380952381 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 27.78, "offset": 0, "label": "rectangular", "acEffectiveFrequency": 248423.92512497256, "effectiveFrequency": 248408.7565140517, "peak": 13.89, "rms": 13.890000000000022, "thd": 0.48331514845248535 }, "harmonics": {"amplitudes": [0.21703125, 17.681745968226167, 0.4340625, 5.884441743008607, 0.4340625, 3.5192857754085134, 0.4340625, 2.50156382659689, 0.4340625, 1.932968090534883, 0.4340625, 1.56850033166751, 0.4340625, 1.313925940874188, 0.4340625, 1.1252647178556892, 0.4340625, 0.9792293094780014, 0.4340625, 0.8623340820515439, 0.4340625, 0.7662274695502621, 0.4340625, 0.6854595927802318, 0.4340625, 0.6163213952979067, 0.4340625, 0.5561996920846202, 0.4340625, 0.5031990666519233, 0.4340625, 0.45591010107099894, 0.4340625, 0.41326185461487136, 0.4340625, 0.3744248874701925, 0.4340625, 0.33874569976854113, 0.4340625, 0.30570130348173624, 0.4340625, 0.2748670467095772, 0.4340625, 0.24589336899764186, 0.4340625, 0.21848870156913297, 0.4340625, 0.1924066733732579, 0.4340625, 0.16743638267205946, 0.4340625, 0.1433948809785266, 0.4340625, 0.12012127132032469, 0.4340625, 0.09747199388796957, 0.4340625, 0.07531698847858959, 0.4340625, 0.05353650312310587, 0.4340625, 0.03201837355771753, 0.4340625, 0.0106556362841701 ], "frequencies": [0, 42000, 84000, 126000, 168000, 210000, 252000, 294000, 336000, 378000, 420000, 462000, 504000, 546000, 588000, 630000, 672000, 714000, 756000, 798000, 840000, 882000, 924000, 966000, 1008000, 1050000, 1092000, 1134000, 1176000, 1218000, 1260000, 1302000, 1344000, 1386000, 1428000, 1470000, 1512000, 1554000, 1596000, 1638000, 1680000, 1722000, 1764000, 1806000, 1848000, 1890000, 1932000, 1974000, 2016000, 2058000, 2100000, 2142000, 2184000, 2226000, 2268000, 2310000, 2352000, 2394000, 2436000, 2478000, 2520000, 2562000, 2604000, 2646000 ] } } }, {"name": "Primary winding excitation", "frequency": 42000, "current": {"waveform": {"ancillaryLabel": null, "data": [0, 0, 0.08, 0, 0 ], "numberPeriods": null, "time": [0, 0.000011904761904761905, 0.000011904761904761905, 0.00002380952380952381, 0.00002380952380952381 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 0.08, "offset": 0, "label": "flybackSecondary", "acEffectiveFrequency": 273605.2255259123, "effectiveFrequency": 239619.35667760245, "peak": 0.07999999999999997, "rms": 0.03304234971971578, "thd": 0.6765065136547783 }, "harmonics": {"amplitudes": [0.020312500000000067, 0.030524888492795912, 0.012737510154434986, 0.008814359426194993, 0.006376435773361427, 0.00522572626420801, 0.004259510523085881, 0.0037274726710446416, 0.0032036443096768643, 0.0029039250980520386, 0.0025722262624417523, 0.002383634614056616, 0.002153058872797908, 0.0020257907808299194, 0.0018552057703927027, 0.0017651543060559615, 0.0016332037060954627, 0.0015673323589111663, 0.0014617999167910994, 0.0014124689465512088, 0.001325847107487929, 0.0012882978141013626, 0.0012157102973274464, 0.001186833188496806, 0.001124970278920515, 0.0011026537403985668, 0.0010491870567694052, 0.0010319511635847236, 0.0009851932793140594, 0.000971973314179004, 0.0009306703387528712, 0.0009206833404708095, 0.0008838834764831795, 0.0008765430049266404, 0.0008435104268187535, 0.0008383704463786768, 0.000808527229199984, 0.0008052441820185946, 0.0007781301537945788, 0.0007764367366278998, 0.0007516811086688038, 0.0007513677793952459, 0.0007286687102576646, 0.0007295704201702419, 0.0007086800435204446, 0.0007106665783623426, 0.000691379870043054, 0.0006943487324152177, 0.0006764951251827419, 0.00068036624010201, 0.0006638032388622281, 0.0006685149910647295, 0.000653123268674609, 0.0006586295307287333, 0.0006443091374032789, 0.0006505770495551611, 0.0006372444738801973, 0.0006442528069273796, 0.0006318386998719282, 0.0006395766820603246, 0.0006280241077351936, 0.0006364906326860852, 0.0006257537477939952, 0.0006349569074192705 ], "frequencies": [0, 42000, 84000, 126000, 168000, 210000, 252000, 294000, 336000, 378000, 420000, 462000, 504000, 546000, 588000, 630000, 672000, 714000, 756000, 798000, 840000, 882000, 924000, 966000, 1008000, 1050000, 1092000, 1134000, 1176000, 1218000, 1260000, 1302000, 1344000, 1386000, 1428000, 1470000, 1512000, 1554000, 1596000, 1638000, 1680000, 1722000, 1764000, 1806000, 1848000, 1890000, 1932000, 1974000, 2016000, 2058000, 2100000, 2142000, 2184000, 2226000, 2268000, 2310000, 2352000, 2394000, 2436000, 2478000, 2520000, 2562000, 2604000, 2646000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-7.14, 7.14, 7.14, -7.14, -7.14 ], "numberPeriods": null, "time": [0, 0, 0.000011904761904761905, 0.000011904761904761905, 0.00002380952380952381 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 14.28, "offset": 0, "label": "rectangular", "acEffectiveFrequency": 248423.9251249721, "effectiveFrequency": 248408.75651405123, "peak": 7.14, "rms": 7.139999999999993, "thd": 0.4833151484524845 }, "harmonics": {"amplitudes": [0.1115625, 9.08910483895859, 0.223125, 3.0248318246998886, 0.223125, 1.8090497074454124, 0.223125, 1.2859010598921377, 0.223125, 0.9936207463224668, 0.223125, 0.8062701488917218, 0.223125, 0.6754090149634054, 0.223125, 0.5784298117703113, 0.223125, 0.5033619344616945, 0.223125, 0.44327324304161464, 0.223125, 0.3938707078897675, 0.223125, 0.3523528792261233, 0.223125, 0.3168131578421202, 0.223125, 0.285908265045658, 0.223125, 0.2586638830737746, 0.223125, 0.23435551631727392, 0.223125, 0.21243265960764457, 0.223125, 0.19246894863478564, 0.223125, 0.1741284590602866, 0.223125, 0.1571423547055144, 0.223125, 0.14129234798462062, 0.223125, 0.1263987512342087, 0.223125, 0.1123116867677183, 0.223125, 0.09890451028690123, 0.223125, 0.08606881009924447, 0.223125, 0.07371054356995504, 0.223125, 0.0617470034000806, 0.223125, 0.05010439426638574, 0.223125, 0.0387158601682599, 0.223125, 0.027519843938011213, 0.223125, 0.016458688783449027, 0.223125, 0.005477411308063118 ], "frequencies": [0, 42000, 84000, 126000, 168000, 210000, 252000, 294000, 336000, 378000, 420000, 462000, 504000, 546000, 588000, 630000, 672000, 714000, 756000, 798000, 840000, 882000, 924000, 966000, 1008000, 1050000, 1092000, 1134000, 1176000, 1218000, 1260000, 1302000, 1344000, 1386000, 1428000, 1470000, 1512000, 1554000, 1596000, 1638000, 1680000, 1722000, 1764000, 1806000, 1848000, 1890000, 1932000, 1974000, 2016000, 2058000, 2100000, 2142000, 2184000, 2226000, 2268000, 2310000, 2352000, 2394000, 2436000, 2478000, 2520000, 2562000, 2604000, 2646000 ] } } } ] } ] }, "magnetic": {"coil": {"functionalDescription": [{"name": "Primary"}, {"name": "Secondary"}, {"name": "Tertiary"} ] } } })";
        json masJson = json::parse(masString);

        if (masJson["inputs"]["designRequirements"]["insulation"]["cti"] == "GROUP_IIIa") {
            masJson["inputs"]["designRequirements"]["insulation"]["cti"] = "GROUP_IIIA";
        }
        if (masJson["inputs"]["designRequirements"]["insulation"]["cti"] == "GROUP_IIIb") {
            masJson["inputs"]["designRequirements"]["insulation"]["cti"] = "GROUP_IIIB";
        }

        OpenMagnetics::Inputs inputs(masJson["inputs"]);
        OpenMagnetics::Mas masMagnetic;
        // inputs.process();

        MagneticAdviser MagneticAdviser;
        // MagneticAdviser.set_interleaving_level(interleavingLevel);
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 1);
        REQUIRE(masMagnetics.size() > 0);

        for (size_t i = 0; i < masMagnetics.size(); ++i){
            auto simulatedMas = masMagnetics[i].first;
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "MagneticSimulator_" + std::to_string(i) + "_" + simulatedMas.get_mutable_magnetic().get_mutable_core().get_shape_name() + ".svg";
            filename = std::filesystem::path(std::regex_replace(std::string(filename), std::regex(" "), "_")).string();
            filename = std::filesystem::path(std::regex_replace(std::string(filename), std::regex("/"), "_")).string();

            outFile.append(filename);
            Painter painter(outFile);
            settings.set_painter_mode(PainterModes::CONTOUR);

            settings.set_painter_number_points_x(20);
            settings.set_painter_number_points_y(20);

            settings.set_painter_include_fringing(false);
            settings.set_painter_mirroring_dimension(0);
            painter.paint_magnetic_field(inputs.get_operating_point(0), simulatedMas.get_mutable_magnetic());
            painter.paint_core(simulatedMas.get_mutable_magnetic());
            painter.paint_bobbin(simulatedMas.get_mutable_magnetic());
            painter.paint_coil_turns(simulatedMas.get_mutable_magnetic());
            painter.export_svg();
        }
    }

    TEST_CASE("MagneticSimulatorJsonLV", "[processor][magnetic-simulator]") {
        json masJson = OpenMagneticsTesting::fixtures::get_json("mas-flyback-3-winding-lv");

        if (masJson["inputs"]["designRequirements"]["insulation"]["cti"] == "GROUP_IIIa") {
            masJson["inputs"]["designRequirements"]["insulation"]["cti"] = "GROUP_IIIA";
        }
        if (masJson["inputs"]["designRequirements"]["insulation"]["cti"] == "GROUP_IIIb") {
            masJson["inputs"]["designRequirements"]["insulation"]["cti"] = "GROUP_IIIB";
        }

        OpenMagnetics::Inputs inputs(masJson["inputs"]);
        OpenMagnetics::Mas masMagnetic;
        // inputs.process();

        MagneticAdviser MagneticAdviser;
        // MagneticAdviser.set_interleaving_level(interleavingLevel);
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 1);
        REQUIRE(masMagnetics.size() > 0);

        for (size_t i = 0; i < masMagnetics.size(); ++i){
            auto simulatedMas = masMagnetics[i].first;
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "MagneticSimulator_" + std::to_string(i) + "_" + simulatedMas.get_mutable_magnetic().get_mutable_core().get_shape_name() + ".svg";
            filename = std::filesystem::path(std::regex_replace(std::string(filename), std::regex(" "), "_")).string();
            filename = std::filesystem::path(std::regex_replace(std::string(filename), std::regex("/"), "_")).string();

            outFile.append(filename);
            Painter painter(outFile);

            settings.set_painter_mode(PainterModes::CONTOUR);
            settings.set_painter_number_points_x(20);
            settings.set_painter_number_points_y(20);
            settings.set_painter_include_fringing(false);
            settings.set_painter_mirroring_dimension(0);
            painter.paint_magnetic_field(inputs.get_operating_point(0), simulatedMas.get_mutable_magnetic());
            painter.paint_core(simulatedMas.get_mutable_magnetic());
            painter.paint_bobbin(simulatedMas.get_mutable_magnetic());
            painter.paint_coil_turns(simulatedMas.get_mutable_magnetic());
            painter.export_svg();
        }
    }

    TEST_CASE("Test_IEEE_Article_Simulation_0", "[processor][magnetic-simulator][heavy]") {

        double temperature = 20;
        std::vector<int64_t> numberTurns = {16, 1};
        std::vector<int64_t> numberParallels = {2, 1};
        std::vector<double> turnsRatios = {16};
        std::string shapeName = "PQ 35/35";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        {
            OpenMagnetics::Wire wire;
            wire.set_strand("Round 0.1 - Grade 1");
            wire.set_nominal_value_outer_diameter(0.001);
            wire.set_number_conductors(45);
            wire.set_type(WireType::LITZ);
            InsulationWireCoating coating;
            coating.set_type(InsulationWireCoatingType::SERVED);
            coating.set_number_layers(1);
            wire.set_coating(coating);
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire;
            wire.set_nominal_value_conducting_width(0.0005);
            wire.set_nominal_value_conducting_height(0.020);
            wire.set_nominal_value_outer_width(0.00053);
            wire.set_nominal_value_outer_height(0.020);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            wire.set_type(WireType::FOIL);
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

        coil.wind({0, 1, 0}, 1);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C90";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.00020);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        json inputsJson;

        inputsJson["operatingPoints"] = json::array();
        json operatingPointJson = json();
        operatingPointJson["name"] = "Nominal";
        operatingPointJson["conditions"] = json();
        operatingPointJson["conditions"]["ambientTemperature"] = 85;

        operatingPointJson["excitationsPerWinding"] = json::array();
        {
            json windingExcitation = json();
            windingExcitation["frequency"] = 115000;
            windingExcitation["current"]["processed"]["label"] = "sinusoidal";
            windingExcitation["current"]["processed"]["peakToPeak"] = 4 * 2;
            windingExcitation["current"]["processed"]["offset"] = 0;
            windingExcitation["current"]["processed"]["dutyCycle"] = 0.5;
            windingExcitation["voltage"]["processed"]["label"] = "rectangular";
            windingExcitation["voltage"]["processed"]["peakToPeak"] = 380;
            windingExcitation["voltage"]["processed"]["offset"] = 0;
            windingExcitation["voltage"]["processed"]["dutyCycle"] = 0.5;
            operatingPointJson["excitationsPerWinding"].push_back(windingExcitation);
        }
        {
            json windingExcitation = json();
            windingExcitation["frequency"] = 115000;
            windingExcitation["current"]["processed"]["label"] = "sinusoidal";
            windingExcitation["current"]["processed"]["peakToPeak"] = 64 * 2;
            windingExcitation["current"]["processed"]["offset"] = 0;
            windingExcitation["current"]["processed"]["dutyCycle"] = 0.5;
            windingExcitation["voltage"]["processed"]["label"] = "rectangular";
            windingExcitation["voltage"]["processed"]["peakToPeak"] = 24;
            windingExcitation["voltage"]["processed"]["offset"] = 0;
            windingExcitation["voltage"]["processed"]["dutyCycle"] = 0.5;
            operatingPointJson["excitationsPerWinding"].push_back(windingExcitation);
        }
        inputsJson["operatingPoints"].push_back(operatingPointJson);

        inputsJson["designRequirements"] = json();
        inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 200e-6;
        inputsJson["designRequirements"]["turnsRatios"] = json::array();
        // A bare number parses into a DimensionWithTolerance with nothing set, and the
        // first resolver call then throws "has neither nominal, minimum nor maximum
        // set" (ABT #190). turnsRatios entries are dimensionWithTolerance objects.
        inputsJson["designRequirements"]["turnsRatios"].push_back(json({{"nominal", 16}}));

        OpenMagnetics::Inputs inputs(inputsJson);



        auto windingLosses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
        json windingLossesJson;
        to_json(windingLossesJson, windingLosses);

        auto simulator = MagneticSimulator();
        auto operatingPoint = inputs.get_operating_point(0);
        auto coreLosses = simulator.calculate_core_losses(operatingPoint, magnetic);
        auto mas = simulator.simulate(inputs, magnetic);
        json coreLossesJson;
        to_json(coreLossesJson, coreLosses);

        check_simulation_outputs_sane(mas);
        CHECK(std::isfinite(windingLosses.get_winding_losses()));
        CHECK(windingLosses.get_winding_losses() > 0);
        CHECK(std::isfinite(coreLosses.get_core_losses()));
        CHECK(coreLosses.get_core_losses() > 0);

        // MagneticAdviser::preview_magnetic(mas);

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_0.mas.json");
            to_file(outFile, mas);
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_0_Waveform_Primary_Current.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);
            painter.paint_waveform(inputs.get_winding_excitation(0, 0).get_current()->get_waveform().value());
            painter.export_svg();
            settings.reset();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_0_Waveform_Primary_Voltage.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);
            painter.paint_waveform(inputs.get_winding_excitation(0, 0).get_voltage()->get_waveform().value());
            painter.export_svg();
            settings.reset();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_0_Waveform_Secondaty_Current.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);
            painter.paint_waveform(inputs.get_winding_excitation(0, 1).get_current()->get_waveform().value());
            painter.export_svg();
            settings.reset();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_0_Waveform_Secondaty_Voltage.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);
            painter.paint_waveform(inputs.get_winding_excitation(0, 1).get_voltage()->get_waveform().value());
            painter.export_svg();
            settings.reset();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_0.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);

            settings.set_painter_mode(PainterModes::CONTOUR);
            settings.set_painter_number_points_x(100);
            settings.set_painter_number_points_y(100);
            // settings.set_painter_maximum_value_colorbar(5000);
            // settings.set_painter_minimum_value_colorbar(500);
            settings.set_painter_logarithmic_scale(true);
            settings.set_painter_include_fringing(true);

            painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            settings.reset();
        }
    }

    TEST_CASE("Test_IEEE_Article_Simulation_1", "[processor][magnetic-simulator][heavy]") {

        double temperature = 20;
        std::vector<int64_t> numberTurns = {18, 11};
        std::vector<int64_t> numberParallels = {1, 2};
        std::vector<double> turnsRatios = {16};
        // Bumped from PQ 65/44 → PQ 65/54: the two thick litz cables (900-strand
        // and 720-strand) no longer fit the 65/44 bobbin window after commit
        // a464a455 clamped bobbin wall/column thicknesses to realistic mm-scale
        // values (previously some configurations got µm-thin walls from the
        // spline interpolator extrapolating outside its training range). PQ 65/54
        // is the smallest standard shape from the probe that lets these wires
        // fit; no waveform/loss assertions in the test depend on the core size.
        std::string shapeName = "PQ 65/54";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        {
            OpenMagnetics::Wire wire;
            wire.set_strand("Round 0.05 - Grade 1");
            wire.set_nominal_value_outer_diameter(OpenMagnetics::Wire::get_outer_diameter_insulated_litz(0.00005, 900, 3, 5.08e-05));
            wire.set_number_conductors(900);
            wire.set_type(WireType::LITZ);
            InsulationWireCoating coating;
            coating.set_type(InsulationWireCoatingType::INSULATED);
            coating.set_number_layers(3);
            coating.set_material("ETFE");
            coating.set_thickness_layers(5.08e-05);
            coating.set_temperature_rating(155);
            coating.set_breakdown_voltage(9000);
            wire.set_coating(coating);
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire;
            wire.set_strand("Round 0.05 - Grade 1");
            wire.set_nominal_value_outer_diameter(OpenMagnetics::Wire::get_outer_diameter_served_litz(0.00005, 720, 1, 1));
            wire.set_number_conductors(800);
            wire.set_type(WireType::LITZ);
            InsulationWireCoating coating;
            coating.set_type(InsulationWireCoatingType::SERVED);
            coating.set_number_layers(1);
            wire.set_coating(coating);
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

        coil.wind({0, 1, 0}, 1);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.00001);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        json inputsJson;

        inputsJson["operatingPoints"] = json::array();
        json operatingPointJson = json();
        operatingPointJson["name"] = "Nominal";
        operatingPointJson["conditions"] = json();
        operatingPointJson["conditions"]["ambientTemperature"] = 20;

        operatingPointJson["excitationsPerWinding"] = json::array();
        // {
        //     json windingExcitation = json();
        //     windingExcitation["frequency"] = 100000;
        //     windingExcitation["current"]["waveform"]["data"] = {0, 17.5, 22.5, -17.5, -22.5, 0};
        //     windingExcitation["current"]["waveform"]["time"] = {0, 0, 0.000005, 0.000005, 0.00001, 0.00001};
        //     windingExcitation["voltage"]["waveform"]["data"] = {800, 800, -800, -800, 800};
        //     windingExcitation["voltage"]["waveform"]["time"] = {0, 0.000005, 0.000005, 0.00001, 0.00001};
        //     operatingPointJson["excitationsPerWinding"].push_back(windingExcitation);
        // }


        {
            json windingExcitation = json();
            windingExcitation["frequency"] = 100000;
            windingExcitation["current"]["processed"]["label"] = "sinusoidal";
            windingExcitation["current"]["processed"]["peakToPeak"] = 22.5 * 2;
            windingExcitation["current"]["processed"]["offset"] = 0;
            windingExcitation["current"]["processed"]["dutyCycle"] = 0.5;
            windingExcitation["voltage"]["processed"]["label"] = "rectangular";
            windingExcitation["voltage"]["processed"]["peakToPeak"] = 1600;
            windingExcitation["voltage"]["processed"]["offset"] = 0;
            windingExcitation["voltage"]["processed"]["dutyCycle"] = 0.5;
            operatingPointJson["excitationsPerWinding"].push_back(windingExcitation);
        }

        inputsJson["operatingPoints"].push_back(operatingPointJson);

        inputsJson["designRequirements"] = json();
        inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 720e-6;
        inputsJson["designRequirements"]["turnsRatios"] = json::array();
        json turnRatioSecondary;
        turnRatioSecondary["minimum"] = 1.6;
        turnRatioSecondary["maximum"] = 1.6;
        turnRatioSecondary["nominal"] = 1.6;
        inputsJson["designRequirements"]["turnsRatios"].push_back(turnRatioSecondary);

        OpenMagnetics::Inputs inputs(inputsJson);

        // auto ea = inputs.get_operating_point(0).get_excitations_per_winding()[0];

        auto windingLosses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
        json windingLossesJson;
        to_json(windingLossesJson, windingLosses);

        auto simulator = MagneticSimulator();
        auto operatingPoint = inputs.get_operating_point(0);
        auto coreLosses = simulator.calculate_core_losses(operatingPoint, magnetic);
        auto mas = simulator.simulate(inputs, magnetic);
        json coreLossesJson;
        to_json(coreLossesJson, coreLosses);

        check_simulation_outputs_sane(mas);
        CHECK(std::isfinite(windingLosses.get_winding_losses()));
        CHECK(windingLosses.get_winding_losses() > 0);
        CHECK(std::isfinite(coreLosses.get_core_losses()));
        CHECK(coreLosses.get_core_losses() > 0);

        // MagneticAdviser::preview_magnetic(mas);

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_1.mas.json");
            to_file(outFile, mas);
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_1_Waveform_Primary_Current.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);
            painter.paint_waveform(inputs.get_winding_excitation(0, 0).get_current()->get_waveform().value());
            painter.export_svg();
            settings.reset();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_1_Waveform_Primary_Voltage.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);
            painter.paint_waveform(inputs.get_winding_excitation(0, 0).get_voltage()->get_waveform().value());
            painter.export_svg();
            settings.reset();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_1_Waveform_Secondaty_Current.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);
            painter.paint_waveform(inputs.get_winding_excitation(0, 1).get_current()->get_waveform().value());
            painter.export_svg();
            settings.reset();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_1_Waveform_Secondaty_Voltage.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);
            painter.paint_waveform(inputs.get_winding_excitation(0, 1).get_voltage()->get_waveform().value());
            painter.export_svg();
            settings.reset();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_1.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);

            settings.set_painter_mode(PainterModes::CONTOUR);
            settings.set_painter_number_points_x(100);
            settings.set_painter_number_points_y(100);
            settings.set_painter_logarithmic_scale(true);
            settings.set_painter_include_fringing(true);

            painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            settings.reset();
        }
    }

    TEST_CASE("Test_IEEE_Article_Simulation_2", "[processor][magnetic-simulator][heavy]") {

        double temperature = 20;
        std::vector<int64_t> numberTurns = {21, 13};
        std::vector<int64_t> numberParallels = {1, 2};
        std::vector<double> turnsRatios = {1.6};
        std::string shapeName = "PQ 60/52";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::INNER_OR_TOP;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        // settings.set_coil_fill_sections_with_margin_tape(true);

        std::vector<OpenMagnetics::Wire> wires;
        {
            OpenMagnetics::Wire wire;
            wire.set_strand("Round 0.071 - Grade 1");
            int numberStrands = 600;
            wire.set_nominal_value_outer_diameter(OpenMagnetics::Wire::get_outer_diameter_insulated_litz(0.000071, numberStrands, 3, 5.08e-05));
            wire.set_number_conductors(numberStrands);
            wire.set_type(WireType::LITZ);
            InsulationWireCoating coating;
            coating.set_type(InsulationWireCoatingType::INSULATED);
            coating.set_number_layers(3);
            coating.set_material("FEP");
            coating.set_thickness_layers(5.08e-05);
            coating.set_temperature_rating(155);
            coating.set_breakdown_voltage(9000);
            wire.set_coating(coating);
            wires.push_back(wire);
        }

        {
            OpenMagnetics::Wire wire;
            wire.set_strand("Round 0.071 - Grade 1");
            int numberStrands = 500;
            wire.set_nominal_value_outer_diameter(OpenMagnetics::Wire::get_outer_diameter_bare_litz(0.000071, numberStrands));
            wire.set_number_conductors(numberStrands);
            wire.set_type(WireType::LITZ);
            InsulationWireCoating coating;
            coating.set_type(InsulationWireCoatingType::BARE);
            wire.set_coating(coating);
            wires.push_back(wire);
        }
        int64_t numberStacks = 1;
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         true,
                                                         numberStacks);

        coil.wind({0, 1}, 2);

        std::string coreMaterial = "3C95F";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.00001);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        json inputsJson;

        inputsJson["operatingPoints"] = json::array();
        json operatingPointJson = json();
        operatingPointJson["name"] = "Nominal";
        operatingPointJson["conditions"] = json();
        operatingPointJson["conditions"]["ambientTemperature"] = 80;

        operatingPointJson["excitationsPerWinding"] = json::array();
        {
            json windingExcitation = json();
            windingExcitation["frequency"] = 100000;
            windingExcitation["current"]["waveform"]["data"] = {0, 17.5, 22.5, -17.5, -22.5, 0};
            windingExcitation["current"]["waveform"]["time"] = {0, 0, 0.000005, 0.000005, 0.00001, 0.00001};
            windingExcitation["voltage"]["waveform"]["data"] = {800, 800, -800, -800, 800};
            windingExcitation["voltage"]["waveform"]["time"] = {0, 0.000005, 0.000005, 0.00001, 0.00001};
            operatingPointJson["excitationsPerWinding"].push_back(windingExcitation);
        }
        inputsJson["operatingPoints"].push_back(operatingPointJson);

        inputsJson["designRequirements"] = json();
        inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 720e-6;
        inputsJson["designRequirements"]["turnsRatios"] = json::array();
        json turnRatioSecondary;
        turnRatioSecondary["minimum"] = 1.6;
        turnRatioSecondary["maximum"] = 1.6;
        turnRatioSecondary["nominal"] = 1.6;
        inputsJson["designRequirements"]["turnsRatios"].push_back(turnRatioSecondary);

        OpenMagnetics::Inputs inputs(inputsJson);



        auto windingLosses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
        json windingLossesJson;
        to_json(windingLossesJson, windingLosses);

        CHECK(std::isfinite(windingLosses.get_winding_losses()));
        CHECK(windingLosses.get_winding_losses() > 0);

        // At 80 C ambient in still air this PQ 60/52 dissipates ~61 W (14 W core, 47 W copper at
        // 80 C), which drives the core to ~270 C; with the core just below 3C95F's Curie
        // temperature (210 C) its losses still drive it to ~300 C. No operating point exists below
        // the ceiling. Rewritten for ABT #1855, Alf approved 2026-10-14: both loss entry points hold
        // the core at the last temperature its material's data covers and give their numbers there,
        // and ResultsAnalyzer reports CORE_ABOVE_MATERIAL_DATA_LIMIT for the operating point.
        auto simulator = MagneticSimulator();
        auto operatingPoint = inputs.get_operating_point(0);
        const auto ceiling = InitialPermeability::get_maximum_ferromagnetic_data_temperature(magnetic.get_core().resolve_material());
        REQUIRE(ceiling);
        CoreLossesOutput heldCoreLosses;
        REQUIRE_NOTHROW(heldCoreLosses = simulator.calculate_core_losses(operatingPoint, magnetic));
        INFO("core losses at " << heldCoreLosses.get_temperature() << " C, ceiling " << ceiling.value() << " C");
        CHECK(heldCoreLosses.get_temperature() <= ceiling.value());
        CHECK(heldCoreLosses.get_temperature() > ceiling.value() - MagneticSimulator::kOperatingTemperatureTolerance);
        OpenMagnetics::Mas simulated;
        REQUIRE_NOTHROW(simulated = simulator.simulate(inputs, magnetic));
        REQUIRE(simulated.get_outputs().size() == 1);
        REQUIRE(simulated.get_outputs()[0].get_core_losses());
        REQUIRE(simulated.get_outputs()[0].get_temperature());
        CHECK(simulated.get_outputs()[0].get_core_losses()->get_temperature() > ceiling.value() - MagneticSimulator::kOperatingTemperatureTolerance);
        CHECK(simulated.get_outputs()[0].get_temperature()->get_maximum_temperature() > ceiling.value());
        size_t dataLimitFindings = 0;
        for (const auto& finding : ResultsAnalyzer::analyze_results(simulated)) {
            INFO(finding.message);
            if (finding.code == ResultFindingCode::CORE_ABOVE_MATERIAL_DATA_LIMIT && finding.operatingPointIndex == 0) {
                ++dataLimitFindings;
            }
        }
        CHECK(dataLimitFindings == 1);


        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_2_Waveform_Primary_Current.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);
            painter.paint_waveform(inputs.get_winding_excitation(0, 0).get_current()->get_waveform().value());
            painter.export_svg();
            settings.reset();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_2_Waveform_Primary_Voltage.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);
            painter.paint_waveform(inputs.get_winding_excitation(0, 0).get_voltage()->get_waveform().value());
            painter.export_svg();
            settings.reset();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_2_Waveform_Secondaty_Current.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);
            painter.paint_waveform(inputs.get_winding_excitation(0, 1).get_current()->get_waveform().value());
            painter.export_svg();
            settings.reset();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_2_Waveform_Secondaty_Voltage.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);
            painter.paint_waveform(inputs.get_winding_excitation(0, 1).get_voltage()->get_waveform().value());
            painter.export_svg();
            settings.reset();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_2.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);

            settings.set_painter_mode(PainterModes::CONTOUR);
            settings.set_painter_number_points_x(100);
            settings.set_painter_number_points_y(100);
            settings.set_painter_logarithmic_scale(true);
            settings.set_painter_include_fringing(true);

            painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            settings.reset();
        }
    }

    // The IEEE article's ETD 29/16/10 3C95 transformer at 100 C ambient. forcedAir: the 2 m/s fan
    // Test_IEEE_Article_Simulation_3 runs it with; without it, the still-air design the fan was added for.
    // `coreMaterial` and `frequency` default to the article's; ABT #1908 drives the same design with other grades
    // and frequencies.
    std::pair<OpenMagnetics::Inputs, OpenMagnetics::Magnetic> build_ieee_article_3(bool forcedAir, std::string coreMaterial = "3C95",
                                                                                   double frequency = 115000) {
        std::vector<int64_t> numberTurns = {32, 2};
        std::vector<int64_t> numberParallels = {1, 1};
        std::vector<double> turnsRatios = {16};
        std::string shapeName = "ETD 29/16/10";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        {
            OpenMagnetics::Wire wire;
            wire.set_strand("Round 0.071 - Grade 1");
            int numberStrands = 110;
            wire.set_nominal_value_outer_diameter(OpenMagnetics::Wire::get_outer_diameter_served_litz(0.000071, numberStrands, 1, 1));
            wire.set_number_conductors(numberStrands);
            wire.set_type(WireType::LITZ);
            InsulationWireCoating coating;
            coating.set_type(InsulationWireCoatingType::SERVED);
            coating.set_number_layers(1);
            wire.set_coating(coating);
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire;
            wire.set_nominal_value_conducting_width(0.0005);
            wire.set_nominal_value_conducting_height(0.018);
            // A foil carries no coating; its insulation is given by insulation layers, so its outer width is its
            // conducting width (Alf, 2026-10-31). Was 0.00053 m on the 0.0005 m sheet.
            wire.set_nominal_value_outer_width(0.0005);
            wire.set_nominal_value_outer_height(0.018);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            wire.set_type(WireType::FOIL);
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

        coil.wind({0, 1}, 1);

        int64_t numberStacks = 1;
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0004);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        json inputsJson;

        inputsJson["operatingPoints"] = json::array();
        json operatingPointJson = json();
        operatingPointJson["name"] = "Nominal";
        operatingPointJson["conditions"] = json();
        operatingPointJson["conditions"]["ambientTemperature"] = 100;
        // Forced air (ABT #838, Alf 2026-10-10; raised from 2 to 4 m/s with Alf's approval of
        // 2026-10-13, "change the fixture so it solves"): this ETD 29 is thermally overloaded at
        // 100 C ambient. main never fed the temperature back into the losses and itself published
        // 224-233 C for it in still air, from the losses at 100 C, without checking it. With the
        // loss-temperature loop it has no operating point in still air (it sheds ~7 W at 3C95's
        // limit and dissipates ~9 W), and at 2 m/s the core is driven to 186.6 C, past the 180 C
        // end of 3C95's measured mu_i(T) table (5.69 W of copper at 100 C under the ABT #1759
        // default winding-loss set). Solved at 100 C ambient: 2.5 m/s still refuses (181.6 C),
        // 3 m/s settles at 174.9 C, 4 m/s at 168.2 C (2.10 W core, 6.10 W copper, hottest point
        // 181.8 C), 5 m/s at 163.3 C. 4 m/s keeps 12 K to the table end, twice what ABT #1759
        // moved this core by; ambient and drive stay the article's.
        if (forcedAir) {
            operatingPointJson["conditions"]["cooling"]["velocity"] = {4.0, 0.0, 0.0};
            operatingPointJson["conditions"]["cooling"]["temperature"] = 100;
        }

        operatingPointJson["excitationsPerWinding"] = json::array();
        {
            json windingExcitation = json();
            windingExcitation["frequency"] = frequency;
            windingExcitation["current"]["processed"]["label"] = "sinusoidal";
            windingExcitation["current"]["processed"]["peakToPeak"] = 4 * 2;
            windingExcitation["current"]["processed"]["offset"] = 0;
            windingExcitation["current"]["processed"]["dutyCycle"] = 0.5;
            windingExcitation["voltage"]["processed"]["label"] = "rectangular";
            windingExcitation["voltage"]["processed"]["peakToPeak"] = 380;
            windingExcitation["voltage"]["processed"]["offset"] = 0;
            windingExcitation["voltage"]["processed"]["dutyCycle"] = 0.5;
            operatingPointJson["excitationsPerWinding"].push_back(windingExcitation);
        }
        {
            json windingExcitation = json();
            windingExcitation["frequency"] = frequency;
            windingExcitation["current"]["processed"]["label"] = "sinusoidal";
            windingExcitation["current"]["processed"]["peakToPeak"] = 64 * 2;
            windingExcitation["current"]["processed"]["offset"] = 0;
            windingExcitation["current"]["processed"]["dutyCycle"] = 0.5;
            windingExcitation["voltage"]["processed"]["label"] = "rectangular";
            windingExcitation["voltage"]["processed"]["peakToPeak"] = 24;
            windingExcitation["voltage"]["processed"]["offset"] = 0;
            windingExcitation["voltage"]["processed"]["dutyCycle"] = 0.5;
            operatingPointJson["excitationsPerWinding"].push_back(windingExcitation);
        }
        inputsJson["operatingPoints"].push_back(operatingPointJson);

        inputsJson["designRequirements"] = json();
        inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 200e-6;
        inputsJson["designRequirements"]["turnsRatios"] = json::array();
        // A bare number parses into a DimensionWithTolerance with nothing set, and the
        // first resolver call then throws "has neither nominal, minimum nor maximum
        // set" (ABT #190). turnsRatios entries are dimensionWithTolerance objects.
        inputsJson["designRequirements"]["turnsRatios"].push_back(json({{"nominal", 16}}));

        return {OpenMagnetics::Inputs(inputsJson), magnetic};
    }

    TEST_CASE("Test_IEEE_Article_Simulation_3", "[processor][magnetic-simulator][heavy]") {

        double temperature = 20;
        auto [inputs, magnetic] = build_ieee_article_3(true);



        auto windingLosses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
        json windingLossesJson;
        to_json(windingLossesJson, windingLosses);

        auto simulator = MagneticSimulator();
        auto operatingPoint = inputs.get_operating_point(0);
        auto coreLosses = simulator.calculate_core_losses(operatingPoint, magnetic);
        auto mas = simulator.simulate(inputs, magnetic);
        json coreLossesJson;
        to_json(coreLossesJson, coreLosses);

        check_simulation_outputs_sane(mas);
        CHECK(std::isfinite(windingLosses.get_winding_losses()));
        CHECK(windingLosses.get_winding_losses() > 0);
        CHECK(std::isfinite(coreLosses.get_core_losses()));
        CHECK(coreLosses.get_core_losses() > 0);

        // MagneticAdviser::preview_magnetic(mas);

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_3.mas.json");
            to_file(outFile, mas);
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_3_Waveform_Primary_Current.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);
            painter.paint_waveform(inputs.get_winding_excitation(0, 0).get_current()->get_waveform().value());
            painter.export_svg();
            settings.reset();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_3_Waveform_Primary_Voltage.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);
            painter.paint_waveform(inputs.get_winding_excitation(0, 0).get_voltage()->get_waveform().value());
            painter.export_svg();
            settings.reset();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_3_Waveform_Secondaty_Current.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);
            painter.paint_waveform(inputs.get_winding_excitation(0, 1).get_current()->get_waveform().value());
            painter.export_svg();
            settings.reset();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_3_Waveform_Secondaty_Voltage.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);
            painter.paint_waveform(inputs.get_winding_excitation(0, 1).get_voltage()->get_waveform().value());
            painter.export_svg();
            settings.reset();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_IEEE_Article_Simulation_3.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);

            settings.set_painter_mode(PainterModes::CONTOUR);
            settings.set_painter_number_points_x(100);
            settings.set_painter_number_points_y(100);
            // settings.set_painter_maximum_value_colorbar(8000);
            // settings.set_painter_minimum_value_colorbar(500);
            settings.set_painter_logarithmic_scale(true);
            settings.set_painter_include_fringing(true);

            painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            settings.reset();
        }
    }

    // ABT #838 x #1739: 3C95's measured mu_i(T) table ends at 180 C, below its 244 C Curie point, and
    // mu_i(T) refuses between the two (MaterialTemperatureOutOfRangeException). In still air at 100 C
    // this ETD 29 has no operating point. The loss-temperature loop used to cap its core iterate only at
    // the Curie point, so the iterate landed in that gap and the MaterialException escaped every
    // adviser's catch (coil adviser's loss ranking included). The loop is capped where mu_i(T) still
    // answers. Rewritten for ABT #1855 (Alf approved 2026-10-14): it no longer refuses there; the core is
    // held at the table's end, 180 C, with its last tabulated permeability, the loss solve and the loss
    // filter give their numbers, and analyze_results reports the data limit for that operating point.
    TEST_CASE("Loss-temperature loop is capped where the core's mu_i(T) data ends, holds the core there and reports the data limit", "[processor][magnetic-simulator]") {
        auto [inputs, magnetic] = build_ieee_article_3(false);
        auto coreMaterial = magnetic.get_mutable_core().resolve_material();
        REQUIRE(coreMaterial.get_name() == "3C95");
        REQUIRE(coreMaterial.get_curie_temperature());
        auto ceiling = InitialPermeability::get_maximum_ferromagnetic_data_temperature(coreMaterial);
        REQUIRE(ceiling);
        // The table ends below Curie: the ceiling is its last point, where mu_i(T) still answers,
        // and just above it mu_i(T) refuses.
        CHECK(ceiling.value() < coreMaterial.get_curie_temperature().value());
        CHECK_NOTHROW(InitialPermeability::get_initial_permeability(coreMaterial, ceiling.value()));
        {
            // The refusal needs the flag off: the default is on since Alf, 2026-10-25.
            const bool allowMaterialDataExtrapolation = settings.get_allow_material_data_extrapolation();
            settings.set_allow_material_data_extrapolation(false);
            CHECK_THROWS_AS(InitialPermeability::get_initial_permeability(coreMaterial, ceiling.value() + 1), MaterialTemperatureOutOfRangeException);
            settings.set_allow_material_data_extrapolation(allowMaterialDataExtrapolation);
        }

        // The loss solve holds the core exactly at the table's end, with the last tabulated permeability
        // (the solved magnetizing current is the real material's there).
        auto operatingPoint = inputs.get_operating_point(0);
        MagneticSimulator simulator;
        CoreLossesOutput coreLosses;
        REQUIRE_NOTHROW(coreLosses = simulator.calculate_core_losses(operatingPoint, magnetic));
        CHECK_THAT(coreLosses.get_temperature(), Catch::Matchers::WithinAbs(ceiling.value(), 1e-9));
        {
            auto operatingPointAtCeiling = inputs.get_operating_point(0);
            operatingPointAtCeiling.get_mutable_conditions().set_ambient_temperature(ceiling.value());
            MagnetizingInductance magnetizingInductanceModel;
            magnetizingInductanceModel.calculate_inductance_and_magnetic_flux_density(magnetic, &operatingPointAtCeiling);
            CHECK_THAT(operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current()->get_processed()->get_peak_to_peak().value(),
                       Catch::Matchers::WithinRel(operatingPointAtCeiling.get_excitations_per_winding()[0].get_magnetizing_current()->get_processed()->get_peak_to_peak().value(), 1e-6));
        }
        // The loss filter the advisers evaluate a candidate with gives its verdict instead of throwing.
        MagneticFilterLosses lossesFilter;
        REQUIRE_NOTHROW(lossesFilter.evaluate_magnetic(&magnetic, &inputs));

        // simulate gives the numbers of the held core, and analyze_results names the data limit.
        OpenMagnetics::Mas mas;
        mas.set_inputs(inputs);
        mas.set_magnetic(magnetic);
        OpenMagnetics::Mas simulated;
        REQUIRE_NOTHROW(simulated = simulator.simulate(mas));
        REQUIRE(simulated.get_outputs().size() == 1);
        CHECK_THAT(simulated.get_outputs()[0].get_core_losses()->get_temperature(), Catch::Matchers::WithinAbs(ceiling.value(), 1e-9));
        auto findings = ResultsAnalyzer::analyze_results(simulated);
        size_t dataLimitFindings = 0;
        for (const auto& finding : findings) {
            INFO(to_string(finding.code) << ": " << finding.message);
            CHECK(finding.operatingPointIndex == 0);
            if (finding.code == ResultFindingCode::CORE_ABOVE_MATERIAL_DATA_LIMIT) {
                ++dataLimitFindings;
                CHECK_THAT(finding.message, Catch::Matchers::ContainsSubstring("last measured permeability"));
                CHECK_THAT(finding.message, Catch::Matchers::ContainsSubstring("the core's real temperature is higher"));
            }
        }
        CHECK(dataLimitFindings == 1);
    }

    // A core held paramagnetic at its Curie point loses almost nothing, so the network the held losses drive comes
    // back COLD: Henry's pick for a 600 W LLC resonant inductor (P 18/11 TM63, 5 turns, 14 uH, 4.17 A rms at 155 kHz,
    // 2026-10-24) loses 23 W at 25 C at 0.37 T, has no fixed point below TM63's 290 C Curie point, and is held there:
    // 0 W of core loss, 1.4 mT and a 33 C hottest point. analyze_results must still say the core left its data,
    // or a design gate reading the hottest point accepts a core sitting at its Curie point; and since ABT #2078 the
    // hottest point itself is published at no less than the Curie point.
    TEST_CASE("A core held paramagnetic at Curie whose network comes back cold is still reported above its data limit", "[processor][magnetic-simulator][abt-1855][curie-cold]") {
        std::ifstream file(OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "curie_held_tm63_p1811_resonant_inductor.json"));
        REQUIRE(file.good());
        auto masJson = nlohmann::json::parse(file);
        OpenMagnetics::Inputs inputs(masJson["inputs"]);
        OpenMagnetics::Magnetic magnetic(masJson["magnetic"]);
        auto coreMaterial = magnetic.get_mutable_core().resolve_material();
        REQUIRE(coreMaterial.get_name() == "TM63");
        auto ceiling = InitialPermeability::get_maximum_ferromagnetic_data_temperature(coreMaterial);
        REQUIRE(ceiling);
        REQUIRE(ceiling.value() == coreMaterial.get_curie_temperature().value());

        MagneticSimulator simulator;
        OpenMagnetics::Mas simulated;
        REQUIRE_NOTHROW(simulated = simulator.simulate(inputs, magnetic));
        REQUIRE(simulated.get_outputs().size() == 1);
        const auto& output = simulated.get_outputs()[0];
        const double coreTemperature = output.get_core_losses()->get_temperature();
        const double hottest = output.get_temperature()->get_maximum_temperature();
        INFO("core losses " << output.get_core_losses()->get_core_losses() << " W at " << coreTemperature << " C, hottest point " << hottest << " C");
        // The premise: held within the tolerance below Curie.
        REQUIRE(coreTemperature > ceiling.value() - MagneticSimulator::kOperatingTemperatureTolerance);
        CHECK(MagneticSimulator::is_core_held_at_data_ceiling(coreMaterial, coreTemperature, hottest));
        // ABT #2078 (Alf, 2026-10-25): the paramagnetic core loses almost nothing and the network comes back cold
        // for it (a 33 C hottest point before), but the published hottest point is never below the material's own
        // limit the core is held at, and the method says it is a floor.
        CHECK(hottest >= ceiling.value());
        CHECK_THAT(output.get_temperature()->get_method_used(), Catch::Matchers::ContainsSubstring("held at its material's data limit"));

        size_t dataLimitFindings = 0;
        for (const auto& finding : ResultsAnalyzer::analyze_results(simulated)) {
            INFO(to_string(finding.code) << ": " << finding.message);
            if (finding.code == ResultFindingCode::CORE_ABOVE_MATERIAL_DATA_LIMIT) {
                ++dataLimitFindings;
                CHECK_THAT(finding.message, Catch::Matchers::ContainsSubstring("no longer magnetic"));
            }
        }
        CHECK(dataLimitFindings == 1);
    }

    // ABT #1855 (Alf, 2026-10-14): the same ETD 29 in still air at 100 C. Its 3C95 mu_i(T) table ends at
    // 180 C, below the 244 C Curie point: the material is still magnetic there. Driven to the table's end,
    // the core is held exactly at its last tabulated temperature and keeps its last tabulated permeability
    // (mu_r = 1 is only for a core above its real Curie temperature).
    TEST_CASE("A core driven to the end of a mu_i(T) table that stops below Curie is held there with its last tabulated permeability", "[processor][magnetic-simulator][abt-1855]") {
        auto [inputs, magnetic] = build_ieee_article_3(false);
        auto coreMaterial = magnetic.get_mutable_core().resolve_material();
        REQUIRE(coreMaterial.get_name() == "3C95");
        auto ceiling = InitialPermeability::get_maximum_ferromagnetic_data_temperature(coreMaterial);
        REQUIRE(ceiling);
        REQUIRE(ceiling.value() < coreMaterial.get_curie_temperature().value());

        auto operatingPoint = inputs.get_operating_point(0);
        MagneticSimulator simulator;
        std::optional<TemperatureOutput> temperature;
        CoreLossesOutput coreLosses;
        REQUIRE_NOTHROW(coreLosses = simulator.calculate_core_losses(operatingPoint, magnetic, nullptr, &temperature));
        CHECK_THAT(coreLosses.get_temperature(), Catch::Matchers::WithinAbs(ceiling.value(), 1e-9));
        REQUIRE(temperature);
        CHECK(temperature->get_maximum_temperature() > ceiling.value());
        CHECK(MagneticSimulator::is_core_held_at_data_ceiling(coreMaterial, coreLosses.get_temperature(), temperature->get_maximum_temperature()));

        // Voltage-driven, the flux density does not depend on the permeability; the magnetizing current
        // the magnetizing model derives (and stamps into the solved excitation) does. The solved one is
        // the real material's at the table's end (its last tabulated permeability), not the paramagnetic
        // core's.
        auto solvedMagnetizingCurrent = operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current();
        REQUIRE(solvedMagnetizingCurrent);
        double solvedPeakToPeak = solvedMagnetizingCurrent->get_processed()->get_peak_to_peak().value();
        auto operatingPointAtCeiling = inputs.get_operating_point(0);
        operatingPointAtCeiling.get_mutable_conditions().set_ambient_temperature(ceiling.value());
        auto paramagneticOperatingPoint = operatingPointAtCeiling;
        MagnetizingInductance magnetizingInductanceModel;
        magnetizingInductanceModel.calculate_inductance_and_magnetic_flux_density(magnetic, &operatingPointAtCeiling);
        magnetizingInductanceModel.calculate_inductance_and_magnetic_flux_density(MagneticSimulator::make_paramagnetic_magnetic(magnetic), &paramagneticOperatingPoint);
        double realPeakToPeak = operatingPointAtCeiling.get_excitations_per_winding()[0].get_magnetizing_current()->get_processed()->get_peak_to_peak().value();
        double paramagneticPeakToPeak = paramagneticOperatingPoint.get_excitations_per_winding()[0].get_magnetizing_current()->get_processed()->get_peak_to_peak().value();
        INFO("solved magnetizing current " << solvedPeakToPeak << " A pk-pk; real material at " << ceiling.value() << " C " << realPeakToPeak
             << " A; mu_r = 1 " << paramagneticPeakToPeak << " A");
        CHECK_THAT(solvedPeakToPeak, Catch::Matchers::WithinRel(realPeakToPeak, 1e-6));
        CHECK(paramagneticPeakToPeak > 1.5 * realPeakToPeak);
    }

    // ABT #1908: the loss <-> temperature solve kept the core below the end of its mu_i(T) / Curie data only,
    // and the other temperature-dependent data of a grade can end lower. DMR71's saturation points stop at
    // 120 C and the monotonic Hermite spline through them, extrapolated, reached zero at 211.2 C, below its
    // 255 C Curie point (its mu_i(T) table runs to 260 C). The IEEE article 3 ETD 29 in still air at 100 C is
    // driven past that: the magnetizing model's magnetisation curve then threw INVALID_INPUT ("no positive
    // saturation field strength"), which ended Henry's grade ranking. Past its data the saturation field
    // strength now goes linearly to zero at Curie (Alf, 2026-10-17), so it sets no ceiling of its own: the core
    // settles inside its data (about 233 C), and analyze_results reports no data limit. Rewritten from the
    // hold at the saturation field strength end (11b51471), Alf approved 2026-10-17.
    TEST_CASE("Loss-temperature loop settles a core past its saturation field strength points inside its data", "[processor][magnetic-simulator][abt-1908]") {
        auto [inputs, magnetic] = build_ieee_article_3(false, "DMR71");
        auto coreMaterial = magnetic.get_mutable_core().resolve_material();
        REQUIRE(coreMaterial.get_name() == "DMR71");
        REQUIRE(coreMaterial.get_curie_temperature());
        const double curieTemperature = coreMaterial.get_curie_temperature().value();
        auto ferromagneticCeiling = InitialPermeability::get_maximum_ferromagnetic_data_temperature(coreMaterial);
        REQUIRE(ferromagneticCeiling);
        REQUIRE(ferromagneticCeiling.value() == curieTemperature);

        auto operatingPoint = inputs.get_operating_point(0);
        MagneticSimulator simulator;
        std::optional<TemperatureOutput> temperature;
        CoreLossesOutput coreLosses;
        REQUIRE_NOTHROW(coreLosses = simulator.calculate_core_losses(operatingPoint, magnetic, nullptr, &temperature));
        const double coreTemperature = coreLosses.get_temperature();
        REQUIRE(temperature);
        INFO("core losses " << coreLosses.get_core_losses() << " W at " << coreTemperature << " C, hottest point "
             << temperature->get_maximum_temperature() << " C");
        // Past where the extrapolated saturation points used to reach zero, and below Curie: inside the data,
        // with a positive saturation field strength, and not held at any ceiling.
        CHECK(coreTemperature > 211.2);
        CHECK(coreTemperature < curieTemperature - MagneticSimulator::kOperatingTemperatureTolerance);
        CHECK(Core::get_magnetic_field_strength_saturation(coreMaterial, coreTemperature) > 0);
        CHECK_NOTHROW(InitialPermeability::get_magnetisation_curve(coreMaterial, coreTemperature, std::nullopt));
        CHECK(!MagneticSimulator::is_core_held_at_data_ceiling(coreMaterial, coreTemperature, temperature->get_maximum_temperature()));

        // The loss filter the advisers rank a candidate with gives its verdict instead of throwing.
        MagneticFilterLosses lossesFilter;
        REQUIRE_NOTHROW(lossesFilter.evaluate_magnetic(&magnetic, &inputs));

        OpenMagnetics::Mas mas;
        mas.set_inputs(inputs);
        mas.set_magnetic(magnetic);
        OpenMagnetics::Mas simulated;
        REQUIRE_NOTHROW(simulated = simulator.simulate(mas));
        REQUIRE(simulated.get_outputs().size() == 1);
        CHECK_THAT(simulated.get_outputs()[0].get_core_losses()->get_temperature(), Catch::Matchers::WithinAbs(coreTemperature, 1e-9));
        std::vector<ResultFinding> findings;
        REQUIRE_NOTHROW(findings = ResultsAnalyzer::analyze_results(simulated));
        for (const auto& finding : findings) {
            INFO(to_string(finding.code) << ": " << finding.message);
            CHECK(finding.operatingPointIndex == 0);
            CHECK(finding.code != ResultFindingCode::CORE_ABOVE_MATERIAL_DATA_LIMIT);
        }
    }

    // ABT #1908: DMR24's Steinmetz range for 400 kHz to 999.3 kHz has ct(T) = 1.137 - 0.0055 T, which stops being
    // positive at 206.9 C, below its 280 C Curie point and before its saturation field strength data end. The
    // same ETD 29 at 500 kHz in still air at 100 C is driven past it: the core-loss model then threw
    // INVALID_LOSS_DATA ("Steinmetz temperature factor ... <= 0"). The core is now held below where ct(T) ends,
    // and analyze_results reports the data limit.
    TEST_CASE("Loss-temperature loop holds the core where its Steinmetz ct(T) data end and reports the data limit", "[processor][magnetic-simulator][abt-1908]") {
        const double frequency = 500000;
        auto [inputs, magnetic] = build_ieee_article_3(false, "DMR24", frequency);
        auto coreMaterial = magnetic.get_mutable_core().resolve_material();
        REQUIRE(coreMaterial.get_name() == "DMR24");
        auto ferromagneticCeiling = InitialPermeability::get_maximum_ferromagnetic_data_temperature(coreMaterial);
        REQUIRE(ferromagneticCeiling);
        auto steinmetzRange = CoreLossesModel::get_steinmetz_coefficients(coreMaterial, frequency);
        REQUIRE(steinmetzRange.get_ct0());
        REQUIRE(steinmetzRange.get_ct1());
        REQUIRE(steinmetzRange.get_ct2().value_or(0) == 0);
        const double temperatureFactorEnd = steinmetzRange.get_ct0().value() / steinmetzRange.get_ct1().value();
        REQUIRE(temperatureFactorEnd < ferromagneticCeiling.value());
        REQUIRE(Core::get_magnetic_field_strength_saturation(coreMaterial, temperatureFactorEnd + MagneticSimulator::kOperatingTemperatureTolerance) > 0);
        // A direct query past the end of ct(T) still refuses.
        try {
            CoreLossesModel::get_temperature_factor(steinmetzRange, temperatureFactorEnd + 1);
            FAIL("ct(T) past the end of its positive range must throw");
        }
        catch (const CalculationException& error) {
            CHECK(error.code() == ErrorCode::INVALID_LOSS_DATA);
        }

        auto operatingPoint = inputs.get_operating_point(0);
        MagneticSimulator simulator;
        std::optional<TemperatureOutput> temperature;
        CoreLossesOutput coreLosses;
        REQUIRE_NOTHROW(coreLosses = simulator.calculate_core_losses(operatingPoint, magnetic, nullptr, &temperature));
        const double heldTemperature = coreLosses.get_temperature();
        INFO("core held at " << heldTemperature << " C, ct(T) ends at " << temperatureFactorEnd << " C");
        CHECK(heldTemperature < temperatureFactorEnd);
        CHECK(heldTemperature > temperatureFactorEnd - MagneticSimulator::kOperatingTemperatureTolerance);
        REQUIRE(temperature);
        CHECK(temperature->get_maximum_temperature() > temperatureFactorEnd);

        MagneticFilterLosses lossesFilter;
        REQUIRE_NOTHROW(lossesFilter.evaluate_magnetic(&magnetic, &inputs));

        OpenMagnetics::Mas mas;
        mas.set_inputs(inputs);
        mas.set_magnetic(magnetic);
        OpenMagnetics::Mas simulated;
        REQUIRE_NOTHROW(simulated = simulator.simulate(mas));
        REQUIRE(simulated.get_outputs().size() == 1);
        std::vector<ResultFinding> findings;
        REQUIRE_NOTHROW(findings = ResultsAnalyzer::analyze_results(simulated));
        size_t dataLimitFindings = 0;
        for (const auto& finding : findings) {
            INFO(to_string(finding.code) << ": " << finding.message);
            if (finding.code == ResultFindingCode::CORE_ABOVE_MATERIAL_DATA_LIMIT) {
                ++dataLimitFindings;
                CHECK(finding.operatingPointIndex == 0);
                CHECK_THAT(finding.message, Catch::Matchers::ContainsSubstring("core-loss data at this frequency is valid for"));
            }
        }
        CHECK(dataLimitFindings == 1);
    }

    // ABT-1908-API-BEGIN
    // ABT #1908: the ceiling's parts are the material's own data ends, found on the models that evaluate them.
    // The saturation field strength part was removed (Alf, 2026-10-17): past its data H_sat goes linearly to zero
    // at Curie, so it never ends below the Curie part.
    TEST_CASE("The core data ceiling covers the Steinmetz ct(T) data ends", "[processor][magnetic-simulator][abt-1908]") {
        {
            SteinmetzCoreLossesMethodRangeDatum linear;
            linear.set_k(1); linear.set_alpha(1); linear.set_beta(2);
            linear.set_ct0(1.137); linear.set_ct1(0.0055); linear.set_ct2(0.0);
            CHECK_THAT(CoreLossesModel::get_temperature_factor_data_end(linear, 25).value(), Catch::Matchers::WithinRel(1.137 / 0.0055, 1e-12));
            SteinmetzCoreLossesMethodRangeDatum parabola;
            parabola.set_k(1); parabola.set_alpha(1); parabola.set_beta(2);
            // ct = (T - 150)(T - 250) / 37500: positive below 150 C, refused between 150 C and 250 C.
            parabola.set_ct0(1.0); parabola.set_ct1(400.0 / 37500); parabola.set_ct2(1.0 / 37500);
            CHECK_THAT(CoreLossesModel::get_temperature_factor_data_end(parabola, 25).value(), Catch::Matchers::WithinRel(150.0, 1e-9));
            CHECK(CoreLossesModel::get_temperature_factor_data_end(parabola, 200).value() == 200);
            CHECK(!CoreLossesModel::get_temperature_factor_data_end(parabola, 260));
            SteinmetzCoreLossesMethodRangeDatum positive;
            positive.set_k(1); positive.set_alpha(1); positive.set_beta(2);
            positive.set_ct0(1.2564); positive.set_ct1(0.0114); positive.set_ct2(4.58e-5);
            CHECK(!CoreLossesModel::get_temperature_factor_data_end(positive, 25));
        }
        {
            // The solve's ceiling names the datum that ends first, on the excitation the solve prices.
            auto ceilingOf = [](std::pair<OpenMagnetics::Inputs, OpenMagnetics::Magnetic> design) {
                auto& [inputs, magnetic] = design;
                auto operatingPoint = inputs.get_operating_point(0);
                MagnetizingInductance magnetizingInductanceModel;
                auto magneticFluxDensity = magnetizingInductanceModel.calculate_inductance_and_magnetic_flux_density(magnetic, &operatingPoint).second;
                auto excitation = operatingPoint.get_excitations_per_winding()[0];
                excitation.set_magnetic_flux_density(magneticFluxDensity);
                CoreLosses coreLossesModel;
                return MagneticSimulator::get_core_data_ceiling(magnetic.get_core(), excitation,
                                                                operatingPoint.get_conditions().get_ambient_temperature(), coreLossesModel);
            };
            auto ceiling = ceilingOf(build_ieee_article_3(false, "DMR24", 500000));
            REQUIRE(ceiling);
            CHECK(ceiling->source == CoreDataCeilingSource::STEINMETZ_TEMPERATURE_FACTOR_END);
            CHECK_THAT(ceiling->temperature, Catch::Matchers::WithinRel(1.1374362023251179 / 0.005497448093004713, 1e-9));
            // At 115 kHz DMR24's range has a ct(T) that stays positive, and the saturation field strength sets no
            // ceiling of its own (it reaches zero at Curie): the ceiling is the Curie temperature.
            auto ceilingLow = ceilingOf(build_ieee_article_3(false, "DMR24"));
            REQUIRE(ceilingLow);
            CHECK(ceilingLow->source == CoreDataCeilingSource::CURIE_TEMPERATURE);
            CHECK(ceilingLow->temperature == find_core_material_by_name("DMR24").get_curie_temperature().value());
            // 3C95 keeps its mu_i(T) table end.
            auto ceiling3C95 = ceilingOf(build_ieee_article_3(false));
            REQUIRE(ceiling3C95);
            CHECK(ceiling3C95->source == CoreDataCeilingSource::PERMEABILITY_TABLE_END);
            CHECK(ceiling3C95->temperature == InitialPermeability::get_maximum_ferromagnetic_data_temperature(find_core_material_by_name("3C95")).value());
        }
    }
    // ABT-1908-API-END

    TEST_CASE("Test_Flyback_Simulation", "[processor][magnetic-simulator]") {

        std::vector<int64_t> numberTurns = {80, 8, 6};
        std::vector<int64_t> numberParallels = {1, 2, 6};
        std::vector<double> turnsRatios = {16, 13};
        std::string shapeName = "ER 28";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round T21A01TXXX-1");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
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

        json inputsJson;

        inputsJson["operatingPoints"] = json::array();
        json operatingPointJson = json();
        operatingPointJson["name"] = "Nominal";
        operatingPointJson["conditions"] = json();
        operatingPointJson["conditions"]["ambientTemperature"] = 100;

        // ABT #190b: the captured waveforms below were scaled by 0.0098 — voltage and
        // current together, so load impedance and L stay consistent. As recorded they
        // put ~1054 V peak across 80 turns of ER 28 at 38 kHz, driving B_peak to
        // 25.5 T (two orders of magnitude past saturation); iGSE then reported 41.9 kW
        // of core losses and tripped the check_simulation_outputs_sane gate. Scaled,
        // B_peak is 0.250 T and core losses 76 mW — an operating point a ferrite
        // flyback can actually have, which is what this test needs in order to
        // exercise the simulator on something meaningful.
        operatingPointJson["excitationsPerWinding"] = json::array();
        {
            json windingExcitation = json();
            windingExcitation["frequency"] = 38000;
            windingExcitation["voltage"]["waveform"]["time"] = {0.0, 3.054077224509666e-10, 9.535015487302436e-10, 1.169532824157234e-9, 1.3855640995842243e-9, 4.696816321468088e-8, 9.319885615605851e-8, 3.410788555624003e-7, 5.070348758027365e-7, 5.956876027992601e-7, 6.870688323048804e-7, 9.126054838506735e-7, 0.0000010532418441536443, 0.0000011444070423838393, 0.000001279858652076569, 0.0000014194148560024064, 0.0000016008811273610902, 0.0000017397892374606484, 0.000001875888940979659, 0.0000019681342955869856, 0.0000020141489572529363, 0.000002151112785873655, 0.0000021966953849887516, 0.000002242710046654702, 0.000002288724708320653, 0.0000023347393699866035, 0.0000027009123818353674, 0.0000027922936113409877, 0.0000028380922417315113, 0.000002884106903397462, 0.0000029299055337879856, 0.0000030670853936841295, 0.000003112884024074653, 0.0000032046973161311275, 0.0000033876758064177968, 0.0000034334744368083204, 0.0000035252877288647947, 0.0000035710863592553183, 0.0000036624675887609386, 0.0000037082662191514623, 0.000003754064849541986, 0.000003845878141598462, 0.000003983058001494606, 0.000004257849783837748, 0.000004303648414228271, 0.000004349447044618795, 0.000004532641566180891, 0.000004578440196571415, 0.000004624238826961939, 0.000004670037457352462, 0.0000048990306093050805, 0.000004944829239695604, 0.0000049906278700861295, 0.000005036426500476653, 0.000005082225130867177, 0.0000051280237612577004, 0.000005173822391648224, 0.000005219621022038748, 0.000005265419652429271, 0.000005311218282819795, 0.000005357016913210319, 0.000005402815543600842, 0.000005448614173991366, 0.00000549441280438189, 0.000005540211434772413, 0.000005586010065162939, 0.0000056318086955534605, 0.000005677607325943986, 0.000005723405956334508, 0.000005769204586725033, 0.000005815003217115555, 0.0000058608018475060805, 0.000005906600477896606, 0.000005952399108287128, 0.000005998197738677653, 0.000006043996369068175, 0.00000622719089063027, 0.000006318788151411317, 0.000006364586781801842, 0.00000645618404258289, 0.000006501982672973415, 0.000006547781303363937, 0.000006593579933754462, 0.000006639594595420413, 0.000006685393225810938, 0.00000673119185620146, 0.000006822573085707079, 0.000006868371716097604, 0.00000696018500815408, 0.000007005983638544602, 0.000007142947467165321, 0.000007233896634120085, 0.000007279479233235185, 0.000007824310109862082, 0.000007869892708977182, 0.000007915475308092275, 0.000008006856537597897, 0.00000819048312171085, 0.000008236497783376797, 0.000008467219185532829, 0.00000865214195729834, 0.000008744387311905666, 0.000008836632666512993, 0.000008882647328178943, 0.000008928661989844894, 0.000009020475281901367, 0.00000915787117307294, 0.000009294835001693654, 0.000009431366767763523, 0.00000947694936687862, 0.00000952231593471829, 0.00000974936480519207, 0.000009794947404307163, 0.000009840313972146837, 0.00000988589657126193, 0.00000993147917037703, 0.000009977061769492127, 0.000010114241629388271, 0.000010160040259778797, 0.000010297652182225791, 0.000010343666843891745, 0.000010389681505557692, 0.000010435696167223646, 0.0000106201868764383, 0.000010712432231045626, 0.000010804677585652956, 0.000010896922940260283, 0.000010942937601926233, 0.000010988952263592184, 0.000011034966925258135, 0.000011126780217314607, 0.000011172578847705132, 0.000011218377478095654, 0.00001126417610848618, 0.000011355557337991802, 0.000011401139937106898, 0.000011446722536221995, 0.000011492305135337092, 0.000011537887734452188, 0.000011583470333567285, 0.000011765368667476821, 0.000011947483032661782, 0.000011993065631776879, 0.000012038648230891979, 0.000012084230830007072, 0.000012175612059512694, 0.00001222141068990322, 0.000012267209320293742, 0.000012313007950684267, 0.000012496634534797216, 0.000012680909212736443, 0.000012726923874402394, 0.000012772938536068345, 0.000012865183890675671, 0.000012911198552341622, 0.000012957213214007576, 0.000013003227875673523, 0.000013095257199005424, 0.000013141271860671375, 0.000013233085152727851, 0.000013278883783118373, 0.000013324682413508898, 0.00001337048104389942, 0.000013416279674289946, 0.000013507660903795564, 0.000013553243502910664, 0.000013644408701140858, 0.000013735573899371051, 0.000013781156498486148, 0.000013826739097601244, 0.000014008853462786206, 0.000014054436061901306, 0.000014100018661016399, 0.000014191399890522021, 0.000014282781120027643, 0.000014465975641589737, 0.000014511774271980263, 0.000014649602225702686, 0.000014695616887368637, 0.000014741631549034587, 0.000014787646210700538, 0.000014833660872366489, 0.00001492569019569839, 0.000014971704857364337, 0.000015017719519030291, 0.000015063734180696238, 0.00001515576350402814, 0.00001520177816569409, 0.000015293591457750565, 0.000015385404749807037, 0.000015431203380197566, 0.000015477002010588088, 0.000015751145699104947, 0.000015796728298220047, 0.000015842310897335147, 0.00001588789349645024, 0.00001593347609556534, 0.000015979058694680434, 0.000016024641293795534, 0.000016070223892910627, 0.000016115806492025727, 0.00001616138909114082, 0.00001620697169025592, 0.000016298352919761542, 0.000016435532779657686, 0.00001648133141004821, 0.000016527130040438737, 0.00001657292867082926, 0.000016710540593276254, 0.000016756555254942204, 0.000016894383208664634, 0.000016986412531996535, 0.000017032427193662486, 0.00001707844185532843, 0.000017124456516994387, 0.00001717047117866033, 0.00001726228447071681, 0.00001730829913238276, 0.000017400112424439234, 0.00001753772434688623, 0.000017583522977276757, 0.0000177207028371729, 0.000017812084066678523, 0.000017857666665793616, 0.00001794904789529924, 0.00001799463049441434, 0.00001804021309352943, 0.00001808579569264453, 0.000018222759521265247, 0.000018268342120380347, 0.00001831392471949544, 0.000018359723349885962, 0.000018405305949001062, 0.000018588284439287728, 0.00001863408306967825, 0.00001867988170006878, 0.00001881749362251578, 0.000018909306914572253, 0.000019047134868294676, 0.000019093149529960627, 0.000019139164191626577, 0.000019276992145349, 0.00001932300680701495, 0.00001941482009907143, 0.000019460834760737374, 0.000019506633391127903, 0.000019552432021518425, 0.000019644245313574898, 0.00001969004394396542, 0.00001973584257435594, 0.000019873022434252092, 0.000019964403663757714, 0.000020009986262872808, 0.00002010136749237843, 0.000020146950091493523, 0.000020283913920114245, 0.000020329496519229338, 0.00002042087774873496, 0.000020512258978240582, 0.000020603640207746204, 0.000020649438838136726, 0.000020695237468527255, 0.000020741036098917777, 0.0000207868347293083, 0.00002083263335969882, 0.000020970245282145822, 0.000021108073235868252, 0.000021337714481647148, 0.00002147554243536957, 0.000021567355727426044, 0.000021659169019482516, 0.000021704967649873038, 0.000021750766280263567, 0.00002179656491065409, 0.00002184236354104461, 0.000021888162171435133, 0.000022025342031331283, 0.000022070924630446377, 0.000022116723260836905, 0.00002220810449034252, 0.00002225368708945762, 0.000022345068318963243, 0.000022482032147583965, 0.00002257341337708958, 0.000022664794606595202, 0.000022710593236985724, 0.000022756391867376246, 0.000022939586388938347, 0.000023021246211049767, 0.00002302146224232519, 0.000023021678273600617, 0.000023021894304876046, 0.000023022110336151468, 0.000023022326367426896, 0.000023022542398702325, 0.000023022758429977754, 0.000023022974461253176, 0.000023023190492528604, 0.000023023406523804033, 0.000023023622555079455, 0.000023023838586354884, 0.000023024054617630312, 0.00002302427064890574, 0.00002306790896654199, 0.00002306812499781742, 0.000023068341029092847, 0.000023068557060368276, 0.000023068773091643698, 0.000023068989122919126, 0.000023069205154194555, 0.000023069421185469977, 0.000023113923628207947, 0.000023114139659483376, 0.000023114355690758804, 0.000023114571722034226, 0.000023159938289873898, 0.000023160154321149326, 0.000023160370352424748, 0.00002320573692026442, 0.00002320595295153985, 0.00002325175158193037, 0.0000232975502123209, 0.00002334334884271142, 0.000023389147473101943, 0.00002343494610349247, 0.000023480744733882994, 0.000023526543364273516, 0.000023572341994664038, 0.000023618140625054566, 0.000023663939255445088, 0.00002370973788583561, 0.000023755536516226132, 0.00002380133514661666, 0.000023847133777007183, 0.000023892932407397705, 0.000023938731037788227, 0.000023984529668178755, 0.000024030328298569277, 0.00002416772418974085, 0.000024305120080912423, 0.000024396717341693467, 0.00002444251597208399, 0.000024488314602474517, 0.00002453411323286504, 0.00002457991186325556, 0.00002462571049364609, 0.000024671509124036612, 0.000024717307754427134, 0.000024763106384817656, 0.000024808905015208184, 0.000024854703645598706, 0.00002490050227598923, 0.00002494630090637975, 0.00002499209953677028, 0.0000250378981671608, 0.000025083696797551323, 0.000025129495427941845, 0.000025175294058332374, 0.000025221092688722896, 0.000025266891319113418, 0.000025312689949503946, 0.000025358488579894468, 0.00002540428721028499, 0.000025450085840675512, 0.00002549588447106604, 0.000026316};
            windingExcitation["voltage"]["waveform"]["data"] = {8.202010613207547, 9.563763757861633, 9.97228970125786, 10.108465015723267, 10.244640330188677, 10.108465015723267, 10.380815644654087, 9.199056603773583, 6.009433962264151, 9.563763757861633, 9.427588443396225, 7.2487834119496855, 9.291413128930817, 9.019062499999997, 7.7934846698113205, 9.155237814465407, 8.065835298742137, 9.019062499999997, 8.202010613207547, 8.610536556603773, 8.74671187106918, 8.474361242138365, 8.338185927672955, 8.474361242138365, 8.610536556603773, 8.74671187106918, 8.610536556603773, 8.474361242138365, 8.474361242138365, 8.610536556603773, 8.610536556603773, 8.474361242138365, 8.474361242138365, 8.610536556603773, 8.474361242138365, 8.474361242138365, 8.610536556603773, 8.610536556603773, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.610536556603773, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.610536556603773, 8.610536556603773, 8.610536556603773, 8.474361242138365, 8.474361242138365, 8.610536556603773, 8.610536556603773, 8.338185927672955, 7.9296599842767295, 7.7934846698113205, 4.797627751572328, 4.661452437106919, 4.525277122641509, 4.389101808176102, 4.661452437106919, 4.797627751572328, 5.887030267295598, 6.976432783018867, 7.384958726415094, 7.7934846698113205, 7.9296599842767295, 8.065835298742137, 8.202010613207547, 8.202010613207547, 7.9296599842767295, 7.384958726415094, 7.2487834119496855, 6.976432783018867, 5.750854952830189, 5.6146796383647795, 5.342329009433963, 5.206153694968554, 5.069978380503144, 4.933803066037735, 4.797627751572328, 4.797627751572328, 4.933803066037735, 5.069978380503144, 5.206153694968554, 5.342329009433963, 6.1593808962264145, 6.5679068396226405, 6.976432783018867, 7.384958726415094, 7.521134040880502, 7.657309355345911, 7.7934846698113205, 7.9296599842767295, 7.9296599842767295, 7.9296599842767295, 7.9296599842767295, 7.7934846698113205, 7.657309355345911, 7.521134040880502, 7.384958726415094, 7.2487834119496855, 7.112608097484276, 6.295556210691824, 5.6146796383647795, 5.4785043238993705, 5.342329009433963, 5.206153694968554, 5.069978380503144, 5.069978380503144, 5.069978380503144, 5.069978380503144, 5.342329009433963, 6.0232055817610055, 6.1593808962264145, 6.295556210691824, 6.7040821540880495, 6.840257468553459, 6.976432783018867, 7.112608097484276, 7.384958726415094, 7.521134040880502, 7.657309355345911, 7.657309355345911, 7.657309355345911, 7.657309355345911, 7.657309355345911, 7.521134040880502, 7.384958726415094, 7.112608097484276, 6.840257468553459, 6.7040821540880495, 6.5679068396226405, 5.887030267295598, 5.750854952830189, 5.6146796383647795, 5.4785043238993705, 5.342329009433963, 5.342329009433963, 5.3423290094339615, 5.614679638364778, 5.750854952830187, 5.8870302672955965, 6.023205581761004, 6.159380896226413, 6.431731525157232, 6.567906839622639, 6.704082154088048, 6.840257468553457, 7.112608097484274, 7.248783411949683, 7.384958726415092, 7.521134040880499, 7.521134040880499, 7.521134040880499, 7.112608097484274, 6.976432783018864, 6.840257468553457, 6.704082154088048, 6.567906839622639, 6.431731525157232, 6.295556210691822, 6.159380896226413, 6.023205581761004, 5.8870302672955965, 5.750854952830187, 5.614679638364778, 5.478504323899368, 5.478504323899368, 5.478504323899368, 5.478504323899368, 5.614679638364778, 5.750854952830187, 6.023205581761004, 6.295556210691822, 6.431731525157232, 6.567906839622639, 6.704082154088048, 6.840257468553457, 6.976432783018864, 7.112608097484274, 7.248783411949683, 7.384958726415092, 7.384958726415092, 7.248783411949683, 7.112608097484274, 6.976432783018864, 6.840257468553457, 6.704082154088048, 6.567906839622639, 6.431731525157232, 6.159380896226413, 6.023205581761004, 5.8870302672955965, 5.8870302672955965, 5.750854952830187, 5.614679638364778, 5.614679638364778, 5.614679638364778, 5.750854952830187, 5.8870302672955965, 6.159380896226413, 6.295556210691822, 6.431731525157232, 6.704082154088048, 6.840257468553457, 6.976432783018864, 7.112608097484274, 7.112608097484274, 7.112608097484274, 7.248783411949683, 7.248783411949683, 7.248783411949683, 7.112608097484274, 6.976432783018864, 6.840257468553457, 6.704082154088048, 6.567906839622639, 6.295556210691822, 6.159380896226413, 6.023205581761004, 5.8870302672955965, 5.750854952830187, 5.750854952830187, 5.750854952830187, 5.750854952830187, 5.750854952830187, 5.750854952830187, 5.8870302672955965, 6.159380896226413, 6.567906839622639, 6.840257468553457, 6.976432783018864, 7.112608097484274, 7.112608097484274, 7.112608097484274, 7.112608097484274, 7.112608097484274, 7.112608097484274, 6.976432783018864, 6.840257468553457, 6.840257468553457, 6.704082154088048, 6.567906839622639, 6.431731525157232, 6.159380896226413, 6.023205581761004, 5.8870302672955965, 5.8870302672955965, 5.8870302672955965, 5.8870302672955965, 0.37703419811320865, 0.24085888364779928, 0.10468356918238997, 0.03149174528301937, 0.16766705974842872, 0.30384237421383803, 0.4400176886792474, 0.5761930031446523, 0.7123683176100616, 0.8485436320754709, 0.9847189465408802, 1.1208942610062895, 1.257069575471699, 1.3932448899371082, 1.5294202044025176, 0.16766705974842872, 0.30384237421383803, 0.4400176886792474, 0.5761930031446523, 0.7123683176100616, 0.8485436320754709, 0.9847189465408802, 1.1208942610062895, 0.30384237421383803, 0.4400176886792474, 0.5761930031446523, 0.7123683176100616, 0.4400176886792474, 0.5761930031446523, 0.7123683176100616, 0.4400176886792474, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523};
            
            // windingExcitation["current"]["waveform"]["data"] = {-0.01, 0.01, -0.01};
            // windingExcitation["current"]["waveform"]["time"] = {0, 0.000026316 / 2, 0.000026316};
            // windingExcitation["current"]["waveform"]["ancillaryLabel"] = "triangular";
            windingExcitation["current"]["waveform"]["time"] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.093188969094081e-8, 3.093188969094081e-8, 3.093188969094081e-8, 3.0931889690942506e-8, 3.0931889690942506e-8, 1.6832778086252022e-7, 1.6832778086252022e-7, 1.6832778086252022e-7, 1.6832778086252022e-7, 2.1412641125304387e-7, 2.599250416435675e-7, 4.0421987784568223e-7, 5.347168239867077e-7, 6.218767351471702e-7, 7.151582454192478e-7, 9.328151027207477e-7, 9.927031278919459e-7, 0.000001212659633536962, 0.0000014303164908384584, 0.000001523598001110531, 0.000001710161021654671, 0.0000018345363686841017, 0.000001927817878956171, 0.000002114380899500316, 0.0000024875069405885978, 0.0000036068850638534516, 0.00002293024708495284, 0.000022976045715343368, 0.000022976045715343368, 0.000022976045715343368, 0.000022976045715343368, 0.00002302184434573389, 0.00002302184434573389, 0.00002302184434573389, 0.00002302184434573389, 0.00002302184434573389, 0.00002302184434573389, 0.00002302184434573389, 0.00002302184434573389, 0.00002302184434573389, 0.00002302184434573389, 0.000023021844345733896, 0.000023021844345733896, 0.000023021844345733896, 0.000023021844345733896, 0.000023021844345733896, 0.000023021844345733896, 0.000023021844345733896, 0.000023021844345733896, 0.000023021844345733896, 0.000023021844345733896, 0.000023021844345733896, 0.000023021844345733896, 0.000023021844345733896, 0.000023021844345733896, 0.000023021844345733896, 0.000023021844345733896, 0.000023021844345733896, 0.000023021844345733896, 0.000023021844345733896, 0.00002302498611882292, 0.00002306764297612441, 0.00002306764297612441, 0.00002306764297612441, 0.00002306764297612441, 0.00002306764297612441, 0.00002306764297612441, 0.00002306764297612441, 0.00002306764297612441, 0.000023067642976124418, 0.000023067642976124418, 0.000023067642976124418, 0.000023067642976124418, 0.000023067642976124418, 0.000023067642976124418, 0.000023067642976124418, 0.000023067642976124418, 0.000023067642976124418, 0.000023067642976124418, 0.000023067642976124418, 0.000023067642976124418, 0.000023067642976124418, 0.000023067642976124418, 0.000023067642976124418, 0.000023067642976124418, 0.000023067642976124418, 0.00002307162687395894, 0.000023113441606514933, 0.000023113441606514933, 0.000023113441606514933, 0.000023113441606514933, 0.000023113441606514933, 0.000023159240236905462, 0.000023159240236905462, 0.000023159240236905462, 0.000023159240236905462, 0.000023159240236905462, 0.000023159240236905462, 0.000023159240236905462, 0.000023205038867295984, 0.000023205038867295984, 0.000023205038867295984, 0.000023205038867295984, 0.000023205038867295984, 0.000023250837497686506, 0.000023250837497686506, 0.000023250837497686506, 0.000023296636128077028, 0.000023296636128077028, 0.000023342434758467557, 0.00002338823338885808, 0.00002338823338885808, 0.0000234340320192486, 0.000023479830649639123, 0.000023708823801591746, 0.000025559133814547507, 0.000025590227651304868, 0.000025632366277993736, 0.000025632366277993736, 0.000026316};
            windingExcitation["current"]["waveform"]["data"] = {0.015202152777778026, 0.031011458333333523, 0.030324097222222422, 0.029636736111111322, 0.02894937500000022, 0.028262013888889118, 0.027574652777777982, 0.026887291666666882, 0.02619993055555578, 0.02551256944444468, 0.024825208333333543, 0.024137847222222442, 0.023450486111111342, 0.02276312500000024, 0.017264236111111362, 0.014514791666666924, 0.013827430555555823, 0.017264236111111362, 0.01657687500000026, 0.015889513888889126, 0.015202152777778026, 0.014514791666666924, 8.020833333362844e-05, -0.0006071527777774725, -0.0012945138888885736, -0.0019818749999996744, -0.0012945138888885736, -0.0012945138888885736, -0.010266666666666317, 0.0048917361111114055, 0.012600000000000269, -0.0013999999999996817, -0.0023333333333330117, 0.0048917361111114055, -0.004199999999999672, 0.0037333333333336334, -0.0023333333333330117, 0.0018666666666669733, -0.0013999999999996817, 0.0018666666666669733, -0.00046666666666635165, 0.0014000000000002734, 1.93493464702042e-17, 8.020833333359363e-05, 0.004891736111111371, 0.0014549305555558304, 0.0007675694444446947, 8.020833333359363e-05, 0.004891736111111371, 0.0042043750000002695, 0.0035170138888891336, 0.0028296527777780327, 0.0021422916666669315, 0.0014549305555558304, 0.0007675694444446947, 8.020833333359363e-05, -0.0006071527777775075, -0.0012945138888886085, -0.0019818749999997096, -0.0026692361111108455, -0.0033565972222219464, -0.004043958333333047, -0.0047313194444441485, -0.0054186805555552844, -0.006106041666666386, -0.006793402777777486, -0.007480763888888587, -0.008168124999999688, -0.008855486111110825, -0.009542847222221925, -0.010230208333333025, -0.010917569444444128, -0.011604930555555263, -0.012292291666666365, -0.012979652777777465, -0.013667013888888566, -0.014354374999999702, -0.030799999999999612, 0.0042043750000002695, 0.0035170138888891336, 0.0028296527777780327, 0.0021422916666669315, 0.0007675694444446947, 8.020833333359363e-05, -0.0006071527777775075, -0.0012945138888886085, -0.0019818749999997096, -0.0026692361111108455, -0.0033565972222219464, -0.004043958333333047, -0.0047313194444441485, -0.0054186805555552844, -0.006106041666666386, -0.006793402777777486, -0.007480763888888587, -0.008168124999999688, -0.008855486111110825, -0.009542847222221925, -0.010230208333333025, -0.010917569444444128, -0.011604930555555263, -0.012292291666666365, -0.012979652777777465, 0.017733333333333514, 0.004891736111111371, 0.0042043750000002695, 0.0007675694444446947, 8.020833333359363e-05, -0.0006071527777775075, 0.004891736111111371, 0.0042043750000002695, 0.0035170138888891336, 0.0028296527777780327, 0.0021422916666669315, 0.0007675694444446947, 8.020833333359363e-05, 0.004891736111111371, 0.0014549305555558304, 0.0007675694444446947, 8.020833333359363e-05, -0.0006071527777775075, 0.0007675694444446947, 8.020833333359363e-05, -0.0006071527777775075, 8.020833333359363e-05, -0.0006071527777775075, 8.020833333359363e-05, 0.0021422916666669315, 0.0014549305555558304, 0.0035170138888891336, 0.0028296527777780327, 0.0042043750000002695, 0.01586666666666689, 0.012133333333333569, 0.027574652777777948, 0.026887291666666847, 0.026199930555555747};
            windingExcitation["current"]["waveform"]["ancillaryLabel"] = "flybackPrimary";
            operatingPointJson["excitationsPerWinding"].push_back(windingExcitation);
        }
        {
            json windingExcitation = json();
            windingExcitation["frequency"] = 38000;
            windingExcitation["voltage"]["waveform"]["time"] = {0.0, 3.054077224509666e-10, 9.535015487302436e-10, 1.169532824157234e-9, 1.3855640995842243e-9, 4.696816321468088e-8, 9.319885615605851e-8, 3.410788555624003e-7, 5.070348758027365e-7, 5.956876027992601e-7, 6.870688323048804e-7, 9.126054838506735e-7, 0.0000010532418441536443, 0.0000011444070423838393, 0.000001279858652076569, 0.0000014194148560024064, 0.0000016008811273610902, 0.0000017397892374606484, 0.000001875888940979659, 0.0000019681342955869856, 0.0000020141489572529363, 0.000002151112785873655, 0.0000021966953849887516, 0.000002242710046654702, 0.000002288724708320653, 0.0000023347393699866035, 0.0000027009123818353674, 0.0000027922936113409877, 0.0000028380922417315113, 0.000002884106903397462, 0.0000029299055337879856, 0.0000030670853936841295, 0.000003112884024074653, 0.0000032046973161311275, 0.0000033876758064177968, 0.0000034334744368083204, 0.0000035252877288647947, 0.0000035710863592553183, 0.0000036624675887609386, 0.0000037082662191514623, 0.000003754064849541986, 0.000003845878141598462, 0.000003983058001494606, 0.000004257849783837748, 0.000004303648414228271, 0.000004349447044618795, 0.000004532641566180891, 0.000004578440196571415, 0.000004624238826961939, 0.000004670037457352462, 0.0000048990306093050805, 0.000004944829239695604, 0.0000049906278700861295, 0.000005036426500476653, 0.000005082225130867177, 0.0000051280237612577004, 0.000005173822391648224, 0.000005219621022038748, 0.000005265419652429271, 0.000005311218282819795, 0.000005357016913210319, 0.000005402815543600842, 0.000005448614173991366, 0.00000549441280438189, 0.000005540211434772413, 0.000005586010065162939, 0.0000056318086955534605, 0.000005677607325943986, 0.000005723405956334508, 0.000005769204586725033, 0.000005815003217115555, 0.0000058608018475060805, 0.000005906600477896606, 0.000005952399108287128, 0.000005998197738677653, 0.000006043996369068175, 0.00000622719089063027, 0.000006318788151411317, 0.000006364586781801842, 0.00000645618404258289, 0.000006501982672973415, 0.000006547781303363937, 0.000006593579933754462, 0.000006639594595420413, 0.000006685393225810938, 0.00000673119185620146, 0.000006822573085707079, 0.000006868371716097604, 0.00000696018500815408, 0.000007005983638544602, 0.000007142947467165321, 0.000007233896634120085, 0.000007279479233235185, 0.000007824310109862082, 0.000007869892708977182, 0.000007915475308092275, 0.000008006856537597897, 0.00000819048312171085, 0.000008236497783376797, 0.000008467219185532829, 0.00000865214195729834, 0.000008744387311905666, 0.000008836632666512993, 0.000008882647328178943, 0.000008928661989844894, 0.000009020475281901367, 0.00000915787117307294, 0.000009294835001693654, 0.000009431366767763523, 0.00000947694936687862, 0.00000952231593471829, 0.00000974936480519207, 0.000009794947404307163, 0.000009840313972146837, 0.00000988589657126193, 0.00000993147917037703, 0.000009977061769492127, 0.000010114241629388271, 0.000010160040259778797, 0.000010297652182225791, 0.000010343666843891745, 0.000010389681505557692, 0.000010435696167223646, 0.0000106201868764383, 0.000010712432231045626, 0.000010804677585652956, 0.000010896922940260283, 0.000010942937601926233, 0.000010988952263592184, 0.000011034966925258135, 0.000011126780217314607, 0.000011172578847705132, 0.000011218377478095654, 0.00001126417610848618, 0.000011355557337991802, 0.000011401139937106898, 0.000011446722536221995, 0.000011492305135337092, 0.000011537887734452188, 0.000011583470333567285, 0.000011765368667476821, 0.000011947483032661782, 0.000011993065631776879, 0.000012038648230891979, 0.000012084230830007072, 0.000012175612059512694, 0.00001222141068990322, 0.000012267209320293742, 0.000012313007950684267, 0.000012496634534797216, 0.000012680909212736443, 0.000012726923874402394, 0.000012772938536068345, 0.000012865183890675671, 0.000012911198552341622, 0.000012957213214007576, 0.000013003227875673523, 0.000013095257199005424, 0.000013141271860671375, 0.000013233085152727851, 0.000013278883783118373, 0.000013324682413508898, 0.00001337048104389942, 0.000013416279674289946, 0.000013507660903795564, 0.000013553243502910664, 0.000013644408701140858, 0.000013735573899371051, 0.000013781156498486148, 0.000013826739097601244, 0.000014008853462786206, 0.000014054436061901306, 0.000014100018661016399, 0.000014191399890522021, 0.000014282781120027643, 0.000014465975641589737, 0.000014511774271980263, 0.000014649602225702686, 0.000014695616887368637, 0.000014741631549034587, 0.000014787646210700538, 0.000014833660872366489, 0.00001492569019569839, 0.000014971704857364337, 0.000015017719519030291, 0.000015063734180696238, 0.00001515576350402814, 0.00001520177816569409, 0.000015293591457750565, 0.000015385404749807037, 0.000015431203380197566, 0.000015477002010588088, 0.000015751145699104947, 0.000015796728298220047, 0.000015842310897335147, 0.00001588789349645024, 0.00001593347609556534, 0.000015979058694680434, 0.000016024641293795534, 0.000016070223892910627, 0.000016115806492025727, 0.00001616138909114082, 0.00001620697169025592, 0.000016298352919761542, 0.000016435532779657686, 0.00001648133141004821, 0.000016527130040438737, 0.00001657292867082926, 0.000016710540593276254, 0.000016756555254942204, 0.000016894383208664634, 0.000016986412531996535, 0.000017032427193662486, 0.00001707844185532843, 0.000017124456516994387, 0.00001717047117866033, 0.00001726228447071681, 0.00001730829913238276, 0.000017400112424439234, 0.00001753772434688623, 0.000017583522977276757, 0.0000177207028371729, 0.000017812084066678523, 0.000017857666665793616, 0.00001794904789529924, 0.00001799463049441434, 0.00001804021309352943, 0.00001808579569264453, 0.000018222759521265247, 0.000018268342120380347, 0.00001831392471949544, 0.000018359723349885962, 0.000018405305949001062, 0.000018588284439287728, 0.00001863408306967825, 0.00001867988170006878, 0.00001881749362251578, 0.000018909306914572253, 0.000019047134868294676, 0.000019093149529960627, 0.000019139164191626577, 0.000019276992145349, 0.00001932300680701495, 0.00001941482009907143, 0.000019460834760737374, 0.000019506633391127903, 0.000019552432021518425, 0.000019644245313574898, 0.00001969004394396542, 0.00001973584257435594, 0.000019873022434252092, 0.000019964403663757714, 0.000020009986262872808, 0.00002010136749237843, 0.000020146950091493523, 0.000020283913920114245, 0.000020329496519229338, 0.00002042087774873496, 0.000020512258978240582, 0.000020603640207746204, 0.000020649438838136726, 0.000020695237468527255, 0.000020741036098917777, 0.0000207868347293083, 0.00002083263335969882, 0.000020970245282145822, 0.000021108073235868252, 0.000021337714481647148, 0.00002147554243536957, 0.000021567355727426044, 0.000021659169019482516, 0.000021704967649873038, 0.000021750766280263567, 0.00002179656491065409, 0.00002184236354104461, 0.000021888162171435133, 0.000022025342031331283, 0.000022070924630446377, 0.000022116723260836905, 0.00002220810449034252, 0.00002225368708945762, 0.000022345068318963243, 0.000022482032147583965, 0.00002257341337708958, 0.000022664794606595202, 0.000022710593236985724, 0.000022756391867376246, 0.000022939586388938347, 0.000023021246211049767, 0.00002302146224232519, 0.000023021678273600617, 0.000023021894304876046, 0.000023022110336151468, 0.000023022326367426896, 0.000023022542398702325, 0.000023022758429977754, 0.000023022974461253176, 0.000023023190492528604, 0.000023023406523804033, 0.000023023622555079455, 0.000023023838586354884, 0.000023024054617630312, 0.00002302427064890574, 0.00002306790896654199, 0.00002306812499781742, 0.000023068341029092847, 0.000023068557060368276, 0.000023068773091643698, 0.000023068989122919126, 0.000023069205154194555, 0.000023069421185469977, 0.000023113923628207947, 0.000023114139659483376, 0.000023114355690758804, 0.000023114571722034226, 0.000023159938289873898, 0.000023160154321149326, 0.000023160370352424748, 0.00002320573692026442, 0.00002320595295153985, 0.00002325175158193037, 0.0000232975502123209, 0.00002334334884271142, 0.000023389147473101943, 0.00002343494610349247, 0.000023480744733882994, 0.000023526543364273516, 0.000023572341994664038, 0.000023618140625054566, 0.000023663939255445088, 0.00002370973788583561, 0.000023755536516226132, 0.00002380133514661666, 0.000023847133777007183, 0.000023892932407397705, 0.000023938731037788227, 0.000023984529668178755, 0.000024030328298569277, 0.00002416772418974085, 0.000024305120080912423, 0.000024396717341693467, 0.00002444251597208399, 0.000024488314602474517, 0.00002453411323286504, 0.00002457991186325556, 0.00002462571049364609, 0.000024671509124036612, 0.000024717307754427134, 0.000024763106384817656, 0.000024808905015208184, 0.000024854703645598706, 0.00002490050227598923, 0.00002494630090637975, 0.00002499209953677028, 0.0000250378981671608, 0.000025083696797551323, 0.000025129495427941845, 0.000025175294058332374, 0.000025221092688722896, 0.000025266891319113418, 0.000025312689949503946, 0.000025358488579894468, 0.00002540428721028499, 0.000025450085840675512, 0.00002549588447106604, 0.000026316};
            windingExcitation["voltage"]["waveform"]["data"] = {8.202010613207547, 9.563763757861633, 9.97228970125786, 10.108465015723267, 10.244640330188677, 10.108465015723267, 10.380815644654087, 9.199056603773583, 6.009433962264151, 9.563763757861633, 9.427588443396225, 7.2487834119496855, 9.291413128930817, 9.019062499999997, 7.7934846698113205, 9.155237814465407, 8.065835298742137, 9.019062499999997, 8.202010613207547, 8.610536556603773, 8.74671187106918, 8.474361242138365, 8.338185927672955, 8.474361242138365, 8.610536556603773, 8.74671187106918, 8.610536556603773, 8.474361242138365, 8.474361242138365, 8.610536556603773, 8.610536556603773, 8.474361242138365, 8.474361242138365, 8.610536556603773, 8.474361242138365, 8.474361242138365, 8.610536556603773, 8.610536556603773, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.610536556603773, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.610536556603773, 8.610536556603773, 8.610536556603773, 8.474361242138365, 8.474361242138365, 8.610536556603773, 8.610536556603773, 8.338185927672955, 7.9296599842767295, 7.7934846698113205, 4.797627751572328, 4.661452437106919, 4.525277122641509, 4.389101808176102, 4.661452437106919, 4.797627751572328, 5.887030267295598, 6.976432783018867, 7.384958726415094, 7.7934846698113205, 7.9296599842767295, 8.065835298742137, 8.202010613207547, 8.202010613207547, 7.9296599842767295, 7.384958726415094, 7.2487834119496855, 6.976432783018867, 5.750854952830189, 5.6146796383647795, 5.342329009433963, 5.206153694968554, 5.069978380503144, 4.933803066037735, 4.797627751572328, 4.797627751572328, 4.933803066037735, 5.069978380503144, 5.206153694968554, 5.342329009433963, 6.1593808962264145, 6.5679068396226405, 6.976432783018867, 7.384958726415094, 7.521134040880502, 7.657309355345911, 7.7934846698113205, 7.9296599842767295, 7.9296599842767295, 7.9296599842767295, 7.9296599842767295, 7.7934846698113205, 7.657309355345911, 7.521134040880502, 7.384958726415094, 7.2487834119496855, 7.112608097484276, 6.295556210691824, 5.6146796383647795, 5.4785043238993705, 5.342329009433963, 5.206153694968554, 5.069978380503144, 5.069978380503144, 5.069978380503144, 5.069978380503144, 5.342329009433963, 6.0232055817610055, 6.1593808962264145, 6.295556210691824, 6.7040821540880495, 6.840257468553459, 6.976432783018867, 7.112608097484276, 7.384958726415094, 7.521134040880502, 7.657309355345911, 7.657309355345911, 7.657309355345911, 7.657309355345911, 7.657309355345911, 7.521134040880502, 7.384958726415094, 7.112608097484276, 6.840257468553459, 6.7040821540880495, 6.5679068396226405, 5.887030267295598, 5.750854952830189, 5.6146796383647795, 5.4785043238993705, 5.342329009433963, 5.342329009433963, 5.3423290094339615, 5.614679638364778, 5.750854952830187, 5.8870302672955965, 6.023205581761004, 6.159380896226413, 6.431731525157232, 6.567906839622639, 6.704082154088048, 6.840257468553457, 7.112608097484274, 7.248783411949683, 7.384958726415092, 7.521134040880499, 7.521134040880499, 7.521134040880499, 7.112608097484274, 6.976432783018864, 6.840257468553457, 6.704082154088048, 6.567906839622639, 6.431731525157232, 6.295556210691822, 6.159380896226413, 6.023205581761004, 5.8870302672955965, 5.750854952830187, 5.614679638364778, 5.478504323899368, 5.478504323899368, 5.478504323899368, 5.478504323899368, 5.614679638364778, 5.750854952830187, 6.023205581761004, 6.295556210691822, 6.431731525157232, 6.567906839622639, 6.704082154088048, 6.840257468553457, 6.976432783018864, 7.112608097484274, 7.248783411949683, 7.384958726415092, 7.384958726415092, 7.248783411949683, 7.112608097484274, 6.976432783018864, 6.840257468553457, 6.704082154088048, 6.567906839622639, 6.431731525157232, 6.159380896226413, 6.023205581761004, 5.8870302672955965, 5.8870302672955965, 5.750854952830187, 5.614679638364778, 5.614679638364778, 5.614679638364778, 5.750854952830187, 5.8870302672955965, 6.159380896226413, 6.295556210691822, 6.431731525157232, 6.704082154088048, 6.840257468553457, 6.976432783018864, 7.112608097484274, 7.112608097484274, 7.112608097484274, 7.248783411949683, 7.248783411949683, 7.248783411949683, 7.112608097484274, 6.976432783018864, 6.840257468553457, 6.704082154088048, 6.567906839622639, 6.295556210691822, 6.159380896226413, 6.023205581761004, 5.8870302672955965, 5.750854952830187, 5.750854952830187, 5.750854952830187, 5.750854952830187, 5.750854952830187, 5.750854952830187, 5.8870302672955965, 6.159380896226413, 6.567906839622639, 6.840257468553457, 6.976432783018864, 7.112608097484274, 7.112608097484274, 7.112608097484274, 7.112608097484274, 7.112608097484274, 7.112608097484274, 6.976432783018864, 6.840257468553457, 6.840257468553457, 6.704082154088048, 6.567906839622639, 6.431731525157232, 6.159380896226413, 6.023205581761004, 5.8870302672955965, 5.8870302672955965, 5.8870302672955965, 5.8870302672955965, 0.37703419811320865, 0.24085888364779928, 0.10468356918238997, 0.03149174528301937, 0.16766705974842872, 0.30384237421383803, 0.4400176886792474, 0.5761930031446523, 0.7123683176100616, 0.8485436320754709, 0.9847189465408802, 1.1208942610062895, 1.257069575471699, 1.3932448899371082, 1.5294202044025176, 0.16766705974842872, 0.30384237421383803, 0.4400176886792474, 0.5761930031446523, 0.7123683176100616, 0.8485436320754709, 0.9847189465408802, 1.1208942610062895, 0.30384237421383803, 0.4400176886792474, 0.5761930031446523, 0.7123683176100616, 0.4400176886792474, 0.5761930031446523, 0.7123683176100616, 0.4400176886792474, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523};

            windingExcitation["current"]["waveform"]["data"] = {-9.8e-05, 9.8e-05, -9.8e-05};
            windingExcitation["current"]["waveform"]["time"] = {0, 0.000026316 / 2, 0.000026316};
            windingExcitation["current"]["waveform"]["ancillaryLabel"] = "triangular";
            // windingExcitation["current"]["waveform"]["time"] = {0.0, 2.4875069405902054e-7, 4.0421987784581436e-7, 5.907828983899358e-7, 8.553072367205054e-7, 0.0000010385017582825898, 0.000001212659633537084, 0.0000013590921710162553, 0.000001523598001110653, 0.0000017101610216547896, 0.0000018034425319268656, 0.000001958911715713646, 0.0000021376688876551624, 0.0000023009439200445712, 0.0000024875069405887164, 0.000002916245604294061, 0.0000047023921895244815, 0.000005022982602258149, 0.00000593895521006862, 0.00000616794836202124, 0.000006305344253192806, 0.000008732671663890568, 0.000008778470294281093, 0.000012488159355913517, 0.000022976045715343462, 0.000023021844345733984, 0.000023021844345733984, 0.000023021844345733984, 0.000023021844345733984, 0.000023021844345733984, 0.000023021844345733984, 0.000023021844345733984, 0.000023067642976124513, 0.000023113441606515035, 0.000023113441606515035, 0.000023113441606515035, 0.000023113441606515035, 0.000023159240236905557, 0.000023159240236905557, 0.000023159240236905557, 0.00002320503886729608, 0.00002329663612807713, 0.00002334243475846765, 0.000023434032019248702, 0.000023846219692763413, 0.000026316};
            // windingExcitation["current"]["waveform"]["data"] = {-2.8671874999996305, -19.722222222221973, -25.277777777777573, -10.83333333333303, -18.823784722221973, -12.686631944444144, -16.3888888888886, -12.686631944444144, -15.277777777777487, -12.777777777777459, -14.166666666666373, -12.222222222221916, -13.09577546296265, -11.666666666666345, -11.944444444444144, -10.231770833333002, -5.731192129629278, -4.91290509259224, -2.4580439814810973, -2.0489004629625924, -1.6397569444440592, -0.003182870369982993, -0.003182870369982993, -0.003182870369982993, 0.40596064814855026, 2.0425347222226264, 1.6333912037041216, 1.2242476851855884, 0.8151041666670835, 0.40596064814855026, -0.003182870369982993, -0.4123263888884878, -1.6397569444440592, 1.2242476851855884, 0.8151041666670835, 0.40596064814855026, -0.4123263888884878, 1.2242476851855884, 0.8151041666670835, 0.40596064814855026, -0.003182870369982993, -0.003182870369982993, -0.003182870369982993, 0.40596064814855026, -0.003182870369982993, 0.27777777777814094};
            // windingExcitation["current"]["waveform"]["ancillaryLabel"] = "flybackSecondary";

            operatingPointJson["excitationsPerWinding"].push_back(windingExcitation);
        }
        {
            json windingExcitation = json();
            windingExcitation["frequency"] = 38000;
            windingExcitation["voltage"]["waveform"]["time"] = {0.0, 3.054077224509666e-10, 9.535015487302436e-10, 1.169532824157234e-9, 1.3855640995842243e-9, 4.696816321468088e-8, 9.319885615605851e-8, 3.410788555624003e-7, 5.070348758027365e-7, 5.956876027992601e-7, 6.870688323048804e-7, 9.126054838506735e-7, 0.0000010532418441536443, 0.0000011444070423838393, 0.000001279858652076569, 0.0000014194148560024064, 0.0000016008811273610902, 0.0000017397892374606484, 0.000001875888940979659, 0.0000019681342955869856, 0.0000020141489572529363, 0.000002151112785873655, 0.0000021966953849887516, 0.000002242710046654702, 0.000002288724708320653, 0.0000023347393699866035, 0.0000027009123818353674, 0.0000027922936113409877, 0.0000028380922417315113, 0.000002884106903397462, 0.0000029299055337879856, 0.0000030670853936841295, 0.000003112884024074653, 0.0000032046973161311275, 0.0000033876758064177968, 0.0000034334744368083204, 0.0000035252877288647947, 0.0000035710863592553183, 0.0000036624675887609386, 0.0000037082662191514623, 0.000003754064849541986, 0.000003845878141598462, 0.000003983058001494606, 0.000004257849783837748, 0.000004303648414228271, 0.000004349447044618795, 0.000004532641566180891, 0.000004578440196571415, 0.000004624238826961939, 0.000004670037457352462, 0.0000048990306093050805, 0.000004944829239695604, 0.0000049906278700861295, 0.000005036426500476653, 0.000005082225130867177, 0.0000051280237612577004, 0.000005173822391648224, 0.000005219621022038748, 0.000005265419652429271, 0.000005311218282819795, 0.000005357016913210319, 0.000005402815543600842, 0.000005448614173991366, 0.00000549441280438189, 0.000005540211434772413, 0.000005586010065162939, 0.0000056318086955534605, 0.000005677607325943986, 0.000005723405956334508, 0.000005769204586725033, 0.000005815003217115555, 0.0000058608018475060805, 0.000005906600477896606, 0.000005952399108287128, 0.000005998197738677653, 0.000006043996369068175, 0.00000622719089063027, 0.000006318788151411317, 0.000006364586781801842, 0.00000645618404258289, 0.000006501982672973415, 0.000006547781303363937, 0.000006593579933754462, 0.000006639594595420413, 0.000006685393225810938, 0.00000673119185620146, 0.000006822573085707079, 0.000006868371716097604, 0.00000696018500815408, 0.000007005983638544602, 0.000007142947467165321, 0.000007233896634120085, 0.000007279479233235185, 0.000007824310109862082, 0.000007869892708977182, 0.000007915475308092275, 0.000008006856537597897, 0.00000819048312171085, 0.000008236497783376797, 0.000008467219185532829, 0.00000865214195729834, 0.000008744387311905666, 0.000008836632666512993, 0.000008882647328178943, 0.000008928661989844894, 0.000009020475281901367, 0.00000915787117307294, 0.000009294835001693654, 0.000009431366767763523, 0.00000947694936687862, 0.00000952231593471829, 0.00000974936480519207, 0.000009794947404307163, 0.000009840313972146837, 0.00000988589657126193, 0.00000993147917037703, 0.000009977061769492127, 0.000010114241629388271, 0.000010160040259778797, 0.000010297652182225791, 0.000010343666843891745, 0.000010389681505557692, 0.000010435696167223646, 0.0000106201868764383, 0.000010712432231045626, 0.000010804677585652956, 0.000010896922940260283, 0.000010942937601926233, 0.000010988952263592184, 0.000011034966925258135, 0.000011126780217314607, 0.000011172578847705132, 0.000011218377478095654, 0.00001126417610848618, 0.000011355557337991802, 0.000011401139937106898, 0.000011446722536221995, 0.000011492305135337092, 0.000011537887734452188, 0.000011583470333567285, 0.000011765368667476821, 0.000011947483032661782, 0.000011993065631776879, 0.000012038648230891979, 0.000012084230830007072, 0.000012175612059512694, 0.00001222141068990322, 0.000012267209320293742, 0.000012313007950684267, 0.000012496634534797216, 0.000012680909212736443, 0.000012726923874402394, 0.000012772938536068345, 0.000012865183890675671, 0.000012911198552341622, 0.000012957213214007576, 0.000013003227875673523, 0.000013095257199005424, 0.000013141271860671375, 0.000013233085152727851, 0.000013278883783118373, 0.000013324682413508898, 0.00001337048104389942, 0.000013416279674289946, 0.000013507660903795564, 0.000013553243502910664, 0.000013644408701140858, 0.000013735573899371051, 0.000013781156498486148, 0.000013826739097601244, 0.000014008853462786206, 0.000014054436061901306, 0.000014100018661016399, 0.000014191399890522021, 0.000014282781120027643, 0.000014465975641589737, 0.000014511774271980263, 0.000014649602225702686, 0.000014695616887368637, 0.000014741631549034587, 0.000014787646210700538, 0.000014833660872366489, 0.00001492569019569839, 0.000014971704857364337, 0.000015017719519030291, 0.000015063734180696238, 0.00001515576350402814, 0.00001520177816569409, 0.000015293591457750565, 0.000015385404749807037, 0.000015431203380197566, 0.000015477002010588088, 0.000015751145699104947, 0.000015796728298220047, 0.000015842310897335147, 0.00001588789349645024, 0.00001593347609556534, 0.000015979058694680434, 0.000016024641293795534, 0.000016070223892910627, 0.000016115806492025727, 0.00001616138909114082, 0.00001620697169025592, 0.000016298352919761542, 0.000016435532779657686, 0.00001648133141004821, 0.000016527130040438737, 0.00001657292867082926, 0.000016710540593276254, 0.000016756555254942204, 0.000016894383208664634, 0.000016986412531996535, 0.000017032427193662486, 0.00001707844185532843, 0.000017124456516994387, 0.00001717047117866033, 0.00001726228447071681, 0.00001730829913238276, 0.000017400112424439234, 0.00001753772434688623, 0.000017583522977276757, 0.0000177207028371729, 0.000017812084066678523, 0.000017857666665793616, 0.00001794904789529924, 0.00001799463049441434, 0.00001804021309352943, 0.00001808579569264453, 0.000018222759521265247, 0.000018268342120380347, 0.00001831392471949544, 0.000018359723349885962, 0.000018405305949001062, 0.000018588284439287728, 0.00001863408306967825, 0.00001867988170006878, 0.00001881749362251578, 0.000018909306914572253, 0.000019047134868294676, 0.000019093149529960627, 0.000019139164191626577, 0.000019276992145349, 0.00001932300680701495, 0.00001941482009907143, 0.000019460834760737374, 0.000019506633391127903, 0.000019552432021518425, 0.000019644245313574898, 0.00001969004394396542, 0.00001973584257435594, 0.000019873022434252092, 0.000019964403663757714, 0.000020009986262872808, 0.00002010136749237843, 0.000020146950091493523, 0.000020283913920114245, 0.000020329496519229338, 0.00002042087774873496, 0.000020512258978240582, 0.000020603640207746204, 0.000020649438838136726, 0.000020695237468527255, 0.000020741036098917777, 0.0000207868347293083, 0.00002083263335969882, 0.000020970245282145822, 0.000021108073235868252, 0.000021337714481647148, 0.00002147554243536957, 0.000021567355727426044, 0.000021659169019482516, 0.000021704967649873038, 0.000021750766280263567, 0.00002179656491065409, 0.00002184236354104461, 0.000021888162171435133, 0.000022025342031331283, 0.000022070924630446377, 0.000022116723260836905, 0.00002220810449034252, 0.00002225368708945762, 0.000022345068318963243, 0.000022482032147583965, 0.00002257341337708958, 0.000022664794606595202, 0.000022710593236985724, 0.000022756391867376246, 0.000022939586388938347, 0.000023021246211049767, 0.00002302146224232519, 0.000023021678273600617, 0.000023021894304876046, 0.000023022110336151468, 0.000023022326367426896, 0.000023022542398702325, 0.000023022758429977754, 0.000023022974461253176, 0.000023023190492528604, 0.000023023406523804033, 0.000023023622555079455, 0.000023023838586354884, 0.000023024054617630312, 0.00002302427064890574, 0.00002306790896654199, 0.00002306812499781742, 0.000023068341029092847, 0.000023068557060368276, 0.000023068773091643698, 0.000023068989122919126, 0.000023069205154194555, 0.000023069421185469977, 0.000023113923628207947, 0.000023114139659483376, 0.000023114355690758804, 0.000023114571722034226, 0.000023159938289873898, 0.000023160154321149326, 0.000023160370352424748, 0.00002320573692026442, 0.00002320595295153985, 0.00002325175158193037, 0.0000232975502123209, 0.00002334334884271142, 0.000023389147473101943, 0.00002343494610349247, 0.000023480744733882994, 0.000023526543364273516, 0.000023572341994664038, 0.000023618140625054566, 0.000023663939255445088, 0.00002370973788583561, 0.000023755536516226132, 0.00002380133514661666, 0.000023847133777007183, 0.000023892932407397705, 0.000023938731037788227, 0.000023984529668178755, 0.000024030328298569277, 0.00002416772418974085, 0.000024305120080912423, 0.000024396717341693467, 0.00002444251597208399, 0.000024488314602474517, 0.00002453411323286504, 0.00002457991186325556, 0.00002462571049364609, 0.000024671509124036612, 0.000024717307754427134, 0.000024763106384817656, 0.000024808905015208184, 0.000024854703645598706, 0.00002490050227598923, 0.00002494630090637975, 0.00002499209953677028, 0.0000250378981671608, 0.000025083696797551323, 0.000025129495427941845, 0.000025175294058332374, 0.000025221092688722896, 0.000025266891319113418, 0.000025312689949503946, 0.000025358488579894468, 0.00002540428721028499, 0.000025450085840675512, 0.00002549588447106604, 0.000026316};
            windingExcitation["voltage"]["waveform"]["data"] = {8.202010613207547, 9.563763757861633, 9.97228970125786, 10.108465015723267, 10.244640330188677, 10.108465015723267, 10.380815644654087, 9.199056603773583, 6.009433962264151, 9.563763757861633, 9.427588443396225, 7.2487834119496855, 9.291413128930817, 9.019062499999997, 7.7934846698113205, 9.155237814465407, 8.065835298742137, 9.019062499999997, 8.202010613207547, 8.610536556603773, 8.74671187106918, 8.474361242138365, 8.338185927672955, 8.474361242138365, 8.610536556603773, 8.74671187106918, 8.610536556603773, 8.474361242138365, 8.474361242138365, 8.610536556603773, 8.610536556603773, 8.474361242138365, 8.474361242138365, 8.610536556603773, 8.474361242138365, 8.474361242138365, 8.610536556603773, 8.610536556603773, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.610536556603773, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.474361242138365, 8.610536556603773, 8.610536556603773, 8.610536556603773, 8.474361242138365, 8.474361242138365, 8.610536556603773, 8.610536556603773, 8.338185927672955, 7.9296599842767295, 7.7934846698113205, 4.797627751572328, 4.661452437106919, 4.525277122641509, 4.389101808176102, 4.661452437106919, 4.797627751572328, 5.887030267295598, 6.976432783018867, 7.384958726415094, 7.7934846698113205, 7.9296599842767295, 8.065835298742137, 8.202010613207547, 8.202010613207547, 7.9296599842767295, 7.384958726415094, 7.2487834119496855, 6.976432783018867, 5.750854952830189, 5.6146796383647795, 5.342329009433963, 5.206153694968554, 5.069978380503144, 4.933803066037735, 4.797627751572328, 4.797627751572328, 4.933803066037735, 5.069978380503144, 5.206153694968554, 5.342329009433963, 6.1593808962264145, 6.5679068396226405, 6.976432783018867, 7.384958726415094, 7.521134040880502, 7.657309355345911, 7.7934846698113205, 7.9296599842767295, 7.9296599842767295, 7.9296599842767295, 7.9296599842767295, 7.7934846698113205, 7.657309355345911, 7.521134040880502, 7.384958726415094, 7.2487834119496855, 7.112608097484276, 6.295556210691824, 5.6146796383647795, 5.4785043238993705, 5.342329009433963, 5.206153694968554, 5.069978380503144, 5.069978380503144, 5.069978380503144, 5.069978380503144, 5.342329009433963, 6.0232055817610055, 6.1593808962264145, 6.295556210691824, 6.7040821540880495, 6.840257468553459, 6.976432783018867, 7.112608097484276, 7.384958726415094, 7.521134040880502, 7.657309355345911, 7.657309355345911, 7.657309355345911, 7.657309355345911, 7.657309355345911, 7.521134040880502, 7.384958726415094, 7.112608097484276, 6.840257468553459, 6.7040821540880495, 6.5679068396226405, 5.887030267295598, 5.750854952830189, 5.6146796383647795, 5.4785043238993705, 5.342329009433963, 5.342329009433963, 5.3423290094339615, 5.614679638364778, 5.750854952830187, 5.8870302672955965, 6.023205581761004, 6.159380896226413, 6.431731525157232, 6.567906839622639, 6.704082154088048, 6.840257468553457, 7.112608097484274, 7.248783411949683, 7.384958726415092, 7.521134040880499, 7.521134040880499, 7.521134040880499, 7.112608097484274, 6.976432783018864, 6.840257468553457, 6.704082154088048, 6.567906839622639, 6.431731525157232, 6.295556210691822, 6.159380896226413, 6.023205581761004, 5.8870302672955965, 5.750854952830187, 5.614679638364778, 5.478504323899368, 5.478504323899368, 5.478504323899368, 5.478504323899368, 5.614679638364778, 5.750854952830187, 6.023205581761004, 6.295556210691822, 6.431731525157232, 6.567906839622639, 6.704082154088048, 6.840257468553457, 6.976432783018864, 7.112608097484274, 7.248783411949683, 7.384958726415092, 7.384958726415092, 7.248783411949683, 7.112608097484274, 6.976432783018864, 6.840257468553457, 6.704082154088048, 6.567906839622639, 6.431731525157232, 6.159380896226413, 6.023205581761004, 5.8870302672955965, 5.8870302672955965, 5.750854952830187, 5.614679638364778, 5.614679638364778, 5.614679638364778, 5.750854952830187, 5.8870302672955965, 6.159380896226413, 6.295556210691822, 6.431731525157232, 6.704082154088048, 6.840257468553457, 6.976432783018864, 7.112608097484274, 7.112608097484274, 7.112608097484274, 7.248783411949683, 7.248783411949683, 7.248783411949683, 7.112608097484274, 6.976432783018864, 6.840257468553457, 6.704082154088048, 6.567906839622639, 6.295556210691822, 6.159380896226413, 6.023205581761004, 5.8870302672955965, 5.750854952830187, 5.750854952830187, 5.750854952830187, 5.750854952830187, 5.750854952830187, 5.750854952830187, 5.8870302672955965, 6.159380896226413, 6.567906839622639, 6.840257468553457, 6.976432783018864, 7.112608097484274, 7.112608097484274, 7.112608097484274, 7.112608097484274, 7.112608097484274, 7.112608097484274, 6.976432783018864, 6.840257468553457, 6.840257468553457, 6.704082154088048, 6.567906839622639, 6.431731525157232, 6.159380896226413, 6.023205581761004, 5.8870302672955965, 5.8870302672955965, 5.8870302672955965, 5.8870302672955965, 0.37703419811320865, 0.24085888364779928, 0.10468356918238997, 0.03149174528301937, 0.16766705974842872, 0.30384237421383803, 0.4400176886792474, 0.5761930031446523, 0.7123683176100616, 0.8485436320754709, 0.9847189465408802, 1.1208942610062895, 1.257069575471699, 1.3932448899371082, 1.5294202044025176, 0.16766705974842872, 0.30384237421383803, 0.4400176886792474, 0.5761930031446523, 0.7123683176100616, 0.8485436320754709, 0.9847189465408802, 1.1208942610062895, 0.30384237421383803, 0.4400176886792474, 0.5761930031446523, 0.7123683176100616, 0.4400176886792474, 0.5761930031446523, 0.7123683176100616, 0.4400176886792474, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523, 0.5761930031446523};

            windingExcitation["current"]["waveform"]["data"] = {-0.00882, 0.00882, -0.00882};
            windingExcitation["current"]["waveform"]["time"] = {0, 0.000026316 / 2, 0.000026316};
            windingExcitation["current"]["waveform"]["ancillaryLabel"] = "triangular";

            operatingPointJson["excitationsPerWinding"].push_back(windingExcitation);
        }
        inputsJson["operatingPoints"].push_back(operatingPointJson);

        inputsJson["designRequirements"] = json();
        inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 1000e-6;
        inputsJson["designRequirements"]["turnsRatios"] = json::array();
        // One ratio per secondary (three windings: 80 : 8 : 6), each a
        // dimensionWithTolerance. A bare number parses into a DimensionWithTolerance
        // with nothing set, and resolve_dimensional_values then throws
        // "has neither nominal, minimum nor maximum set" (ABT #190).
        inputsJson["designRequirements"]["turnsRatios"].push_back(json({{"nominal", 80.0 / 8}}));
        inputsJson["designRequirements"]["turnsRatios"].push_back(json({{"nominal", 80.0 / 6}}));

        OpenMagnetics::Inputs inputs(inputsJson);

        inputs.process();

        auto simulator = MagneticSimulator();
        auto operatingPoint = inputs.get_operating_point(0);
        // auto coreLosses = simulator.calculate_core_losses(operatingPoint, magnetic);
        auto mas = simulator.simulate(inputs, magnetic);
        // json coreLossesJson;
        // to_json(coreLossesJson, coreLosses);

        check_simulation_outputs_sane(mas);

        // MagneticAdviser::preview_magnetic(mas);

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Flyback_Simulation.mas.json");
            to_file(outFile, mas);
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_Flyback_Simulation_Waveform_Primary_Current.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);
            painter.paint_waveform(inputs.get_winding_excitation(0, 0).get_current()->get_waveform().value());
            painter.export_svg();
            settings.reset();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Flyback_Simulation_Waveform_Primary_Voltage.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);
            painter.paint_waveform(inputs.get_winding_excitation(0, 0).get_voltage()->get_waveform().value());
            painter.export_svg();
            settings.reset();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Flyback_Simulation_Waveform_Secondaty_Current.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);
            painter.paint_waveform(inputs.get_winding_excitation(0, 1).get_current()->get_waveform().value());
            painter.export_svg();
            settings.reset();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Flyback_Simulation_Waveform_Secondaty_Voltage.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);
            painter.paint_waveform(inputs.get_winding_excitation(0, 1).get_voltage()->get_waveform().value());
            painter.export_svg();
            settings.reset();
        }

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Flyback_Simulation.svg");
            std::filesystem::remove(outFile); 
            Painter painter(outFile);

            settings.set_painter_mode(PainterModes::CONTOUR);
            settings.set_painter_number_points_x(100);
            settings.set_painter_number_points_y(100);
            settings.set_painter_maximum_value_colorbar(8000);
            settings.set_painter_minimum_value_colorbar(500);
            settings.set_painter_logarithmic_scale(false);
            settings.set_painter_include_fringing(false);

            painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            settings.reset();
        }
    }

    TEST_CASE("Test_Simulator_Web_0", "[processor][magnetic-simulator][bug]") {

        auto json_path_1084 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_simulator_web_0_1084.json");
        std::ifstream json_file_1084(json_path_1084);
        OpenMagnetics::Magnetic magnetic(json::parse(json_file_1084));
        OpenMagnetics::Inputs inputs(json::parse(R"({"designRequirements":{"magnetizingInductance":{"nominal":0.0001},"name":"My Design Requirements","turnsRatios":[{"nominal":1}]},"operatingPoints":[{"name":"Operating Point No. 1","conditions":{"ambientTemperature":42},"excitationsPerWinding":[{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"triangular","acEffectiveFrequency":110746.40291779551,"effectiveFrequency":110746.40291779551,"peak":5,"rms":2.8874560332150576,"thd":0.12151487440704967},"harmonics":{"amplitudes":[1.1608769501236793e-14,4.05366124583194,1.787369544444173e-15,0.4511310569983995,9.749053004706756e-16,0.16293015292554872,4.036157626725542e-16,0.08352979924600704,3.4998295008010614e-16,0.0508569581336163,3.1489164048780735e-16,0.034320410449418075,3.142469873118059e-16,0.024811988673843106,2.3653352035940994e-16,0.018849001010678823,2.9306524147249266e-16,0.014866633059596499,1.796485796132569e-16,0.012077180559557785,1.6247782523152451e-16,0.010049063750920609,1.5324769149805092e-16,0.008529750975091871,1.0558579733068502e-16,0.007363501410705499,7.513269775674661e-17,0.006450045785294609,5.871414177162291e-17,0.005722473794997712,9.294731722001391e-17,0.005134763398167541,1.194820309200107e-16,0.004654430423785411,8.2422739080512e-17,0.004258029771397032,9.5067306351894e-17,0.0039283108282380024,1.7540347128474968e-16,0.0036523670873925395,9.623794010508822e-17,0.0034204021424253787,1.4083470894369491e-16,0.0032248884817922927,1.4749333016985644e-16,0.0030599828465501895,1.0448590642474364e-16,0.002921112944200383,7.575487373767413e-18,0.002804680975178716,7.419510610361002e-17,0.0027078483284668376,3.924741709073613e-17,0.0026283777262804953,2.2684279102637236e-17,0.0025645167846443107,8.997077625295079e-17,0.0025149120164513483,7.131074184849219e-17,0.0024785457043284276,9.354417496250849e-17,0.0024546904085875065,1.2488589642405877e-16,0.0024428775264784264],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}},"voltage":{"waveform":{"data":[-20.5,70.5,70.5,-20.5,-20.5],"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"rectangular","acEffectiveFrequency":591485.5360118389,"effectiveFrequency":553357.3374711702,"peak":70.5,"rms":51.572309672924284,"thd":0.4833151484524849},"harmonics":{"amplitudes":[24.2890625,57.92076613061847,1.421875,19.27588907896988,1.421875,11.528257939603122,1.421875,8.194467538528329,1.421875,6.331896912839248,1.421875,5.137996046859012,1.421875,4.304077056139349,1.421875,3.6860723299088454,1.421875,3.207698601961777,1.421875,2.8247804703632298,1.421875,2.509960393415185,1.421875,2.2453859950684323,1.421875,2.01890737840567,1.421875,1.8219644341144872,1.421875,1.6483482744897402,1.421875,1.4934420157473332,1.421875,1.3537375367153817,1.421875,1.2265178099275544,1.421875,1.1096421410704556,1.421875,1.0013973584174929,1.421875,0.9003924136274832,1.421875,0.8054822382572133,1.421875,0.7157117294021269,1.421875,0.6302738400635857,1.421875,0.5484777114167545,1.421875,0.46972405216147894,1.421875,0.3934858059809043,1.421875,0.31929270856030145,1.421875,0.24671871675852053,1.421875,0.17537155450693565,1.421875,0.10488380107099537,1.421875,0.034905072061178544],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}}},{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"ancillaryLabel":null,"data":[-5,5,-5],"numberPeriods":null,"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"triangular","acEffectiveFrequency":110746.40291779551,"effectiveFrequency":110746.40291779551,"peak":5,"rms":2.8874560332150576,"thd":0.12151487440704967},"harmonics":{"amplitudes":[1.1608769501236793e-14,4.05366124583194,1.787369544444173e-15,0.4511310569983995,9.749053004706756e-16,0.16293015292554872,4.036157626725542e-16,0.08352979924600704,3.4998295008010614e-16,0.0508569581336163,3.1489164048780735e-16,0.034320410449418075,3.142469873118059e-16,0.024811988673843106,2.3653352035940994e-16,0.018849001010678823,2.9306524147249266e-16,0.014866633059596499,1.796485796132569e-16,0.012077180559557785,1.6247782523152451e-16,0.010049063750920609,1.5324769149805092e-16,0.008529750975091871,1.0558579733068502e-16,0.007363501410705499,7.513269775674661e-17,0.006450045785294609,5.871414177162291e-17,0.005722473794997712,9.294731722001391e-17,0.005134763398167541,1.194820309200107e-16,0.004654430423785411,8.2422739080512e-17,0.004258029771397032,9.5067306351894e-17,0.0039283108282380024,1.7540347128474968e-16,0.0036523670873925395,9.623794010508822e-17,0.0034204021424253787,1.4083470894369491e-16,0.0032248884817922927,1.4749333016985644e-16,0.0030599828465501895,1.0448590642474364e-16,0.002921112944200383,7.575487373767413e-18,0.002804680975178716,7.419510610361002e-17,0.0027078483284668376,3.924741709073613e-17,0.0026283777262804953,2.2684279102637236e-17,0.0025645167846443107,8.997077625295079e-17,0.0025149120164513483,7.131074184849219e-17,0.0024785457043284276,9.354417496250849e-17,0.0024546904085875065,1.2488589642405877e-16,0.0024428775264784264],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}},"voltage":{"waveform":{"ancillaryLabel":null,"data":[-50,50,50,-50,-50],"numberPeriods":null,"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"rectangular","acEffectiveFrequency":591485.536011839,"effectiveFrequency":591449.4202715514,"peak":50,"rms":50,"thd":0.48331514845248497},"harmonics":{"amplitudes":[0.78125,63.64919355013018,1.5625,21.18229569117569,1.5625,12.668415318245188,1.5625,9.004909382998164,1.5625,6.958128475647527,1.5625,5.646149502042871,1.5625,4.729755006746538,1.5625,4.050628933965765,1.5625,3.524943518639316,1.5625,3.104154363036517,1.5625,2.7581982345221827,1.5625,2.467457137437843,1.5625,2.2185795367095267,1.5625,2.0021587188071255,1.5625,1.8113717302085082,1.5625,1.6411450722498175,1.5625,1.487623666720196,1.5625,1.3478217691511587,1.5625,1.2193869682092893,1.5625,1.100436657601639,1.5625,0.9894422127774558,1.5625,0.8851453167661671,1.5625,0.7864964059364037,1.5625,0.6926086154544899,1.5625,0.60272275979863,1.5625,0.5161802771005264,1.5625,0.43240198459440116,1.5625,0.3508711083080249,1.5625,0.27111946896540395,1.5625,0.192715993963664,1.5625,0.11525692425384548,1.5625,0.03835722204524927],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}}}]}]})"));

        auto simulator = MagneticSimulator();
        auto mas = simulator.simulate(inputs, magnetic);

        check_simulation_outputs_sane(mas);

        // MagneticAdviser::preview_magnetic(mas);
    }

    TEST_CASE("Test_Simulator_Web_1", "[processor][magnetic-simulator][bug]") {
    // ABT #840: the captured payload carried MORE excitations than the magnetic has windings —
    // op0 held 12 (resp. 6) for a 2-winding magnetic, the surplus being duplicate,
    // under-specified copies (waveform + a partial processed block, no peak/rms/harmonics),
    // and in the web_1 case at the WRONG frequency for an operating point named "25kHz".
    // MKF never read them: every consumer indexes excitationsPerWinding[windingIndex] for
    // windingIndex < numberWindings, so the surplus was dead weight the old code silently
    // ignored. The 0d388fa6 guard now refuses it, correctly — an operating point that supplies
    // an ambiguous number of excitations cannot be simulated without guessing which ones count.
    // The fixtures are trimmed to the excitations MKF actually consumed, so these repros are
    // byte-for-byte the same simulation they always were; only the ignored padding is gone.

        auto json_path_1095 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_simulator_web_1_1095.json");
        std::ifstream json_file_1095(json_path_1095);
        OpenMagnetics::Magnetic magnetic(json::parse(json_file_1095));
        auto json_path_1096 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_simulator_web_1_1096.json");
        std::ifstream json_file_1096(json_path_1096);
        OpenMagnetics::Inputs inputs(json::parse(json_file_1096));

        auto simulator = MagneticSimulator();
        auto mas = simulator.simulate(inputs, magnetic);

        check_simulation_outputs_sane(mas);

        // MagneticAdviser::preview_magnetic(mas);
    }

    TEST_CASE("Test_Simulator_Web_2", "[processor][magnetic-simulator][bug]") {
    // ABT #840: the captured payload carried MORE excitations than the magnetic has windings —
    // op0 held 12 (resp. 6) for a 2-winding magnetic, the surplus being duplicate,
    // under-specified copies (waveform + a partial processed block, no peak/rms/harmonics),
    // and in the web_1 case at the WRONG frequency for an operating point named "25kHz".
    // MKF never read them: every consumer indexes excitationsPerWinding[windingIndex] for
    // windingIndex < numberWindings, so the surplus was dead weight the old code silently
    // ignored. The 0d388fa6 guard now refuses it, correctly — an operating point that supplies
    // an ambiguous number of excitations cannot be simulated without guessing which ones count.
    // The fixtures are trimmed to the excitations MKF actually consumed, so these repros are
    // byte-for-byte the same simulation they always were; only the ignored padding is gone.

        // ABT #1567: the core material is DMR95 (Steinmetz fitted 25 kHz to 500 kHz), not the capture's
        // DMR51W, whose data start at 500 kHz: at this fixture's 100 kHz its losses would be an
        // extrapolation, which MKF refuses since ABT #1456. The contiguous-winding repro is unchanged.
        auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "bug_contiguous.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        auto simulator = MagneticSimulator();
        mas = simulator.simulate(mas.get_inputs(), mas.get_magnetic());

        check_simulation_outputs_sane(mas);

        // MagneticAdviser::preview_magnetic(mas);
    }

    TEST_CASE("Test_Simulator_Web_3", "[processor][magnetic-simulator][bug]") {
        SKIP("Test missing file");
        auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "bug_partialWinding_missing.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        auto simulator = MagneticSimulator();
        mas = simulator.simulate(mas.get_inputs(), mas.get_magnetic());

        check_simulation_outputs_sane(mas);

        // MagneticAdviser::preview_magnetic(mas);
    }

    // =========================================================================
    // Example MAS integration tests
    // Each test loads a MAS JSON from examples/, autocompletes it, simulates,
    // and verifies that core losses, winding losses, temperature, and
    // inductance are present and physically reasonable.
    // =========================================================================

    static std::filesystem::path get_examples_dir() {
        return std::filesystem::path{std::source_location::current().file_name()}
            .parent_path().append("..").append("MAS").append("examples");
    }

    static void run_example_simulation(const std::string& filename) {
        auto path = get_examples_dir() / filename;
        REQUIRE(std::filesystem::exists(path));

        auto mas = OpenMagneticsTesting::mas_loader(path.string());

        auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
        auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);

        REQUIRE(magnetic.get_core().get_processed_description());

        mas.set_magnetic(magnetic);
        mas.set_inputs(inputs);

        MagneticSimulator simulator;
        mas = simulator.simulate(mas.get_inputs(), mas.get_magnetic());

        // Verify outputs exist
        REQUIRE(mas.get_outputs().size() > 0);

        for (size_t opIdx = 0; opIdx < mas.get_outputs().size(); ++opIdx) {
            auto& output = mas.get_outputs()[opIdx];

            // Inductance must be present
            REQUIRE(output.get_inductance());
            auto magInd = output.get_inductance()->get_magnetizing_inductance().get_magnetizing_inductance();
            double inductanceValue = resolve_dimensional_values(magInd);
            REQUIRE(inductanceValue > 0);

            // Core losses must be present and positive
            REQUIRE(output.get_core_losses());
            double coreLosses = output.get_core_losses()->get_core_losses();
            REQUIRE(coreLosses >= 0);

            // Winding losses must be present and positive
            REQUIRE(output.get_winding_losses());
            double windingLosses = output.get_winding_losses()->get_winding_losses();
            REQUIRE(windingLosses >= 0);

            // The core-loss iteration must return a usable core temperature: finite, and above
            // the ambient it started from, since a powered core cannot be colder than the room.
            //
            // This was a flat `< 200 C` ceiling. It passed only because the temperature fed to the
            // core-loss model ignored the winding's heat (ABT #838): counting it puts
            // 09_planar_xfmr_er2510_3c94 at 213.2 C. That design's 188 K rise is dominated not by
            // the thermal network -- 1.56x OMFEM's 2D FEM at identical losses, inside the FEM
            // battery's documented spread -- but by MKF's own losses running 4.4x (core) and 2.1x
            // (copper) above the FEM's for this example, which is ABT #1262. A fixed ceiling here
            // measures that loss discrepancy rather than anything about the simulator, and would
            // have to move again every time the thermal or loss models get more honest.
            const double coreTemperature = output.get_core_losses()->get_temperature();
            const double ambientTemperature =
                mas.get_inputs().get_operating_points()[opIdx].get_conditions().get_ambient_temperature();
            REQUIRE(std::isfinite(coreTemperature));
            REQUIRE(coreTemperature >= ambientTemperature);
        }
    }

    TEST_CASE("Example_01_simple_inductor", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("01_simple_inductor_etd34_n87.json");
    }
    TEST_CASE("Example_02_flyback", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("02_flyback_efd25_3c95.json");
    }
    TEST_CASE("Example_03_buck_inductor", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("03_buck_inductor_pq3230_n95.json");
    }
    TEST_CASE("Example_04_forward_xfmr", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("04_forward_xfmr_e3216_n87.json");
    }
    TEST_CASE("Example_05_pfc_inductor", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("05_pfc_inductor_t4020_hf60.json");
    }
    TEST_CASE("Example_06_llc_xfmr", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("06_llc_xfmr_eq4128_3c97.json");
    }
    TEST_CASE("Example_07_cmc", "[example][magnetic-simulator][smoke-test]") {
        // ABT #1848: the example states a proper common-mode choke operating point at 150 kHz, inside the fitted span of its
        // ferrite's loss data, so it runs like the other examples. The refusal of a frequency outside a material's fitted span
        // (ABT #1456) stays guarded in TestCoreLosses.cpp and TestMaterialDataExtrapolation.cpp.
        run_example_simulation("07_cmc_t2515_w800.json");
    }
    TEST_CASE("Example_08_pq_inductor", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("08_pq_inductor_pq5050_n27.json");
    }
    TEST_CASE("Example_09_planar_xfmr", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("09_planar_xfmr_er2510_3c94.json");
    }
    TEST_CASE("Example_10_emi_filter", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("10_emi_filter_ep13_3f46.json");
    }
    TEST_CASE("Example_11_pushpull", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("11_pushpull_etd49_tp4a.json");
    }
    TEST_CASE("Example_12_boost_inductor", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("12_boost_inductor_t5026_26.json");
    }
    TEST_CASE("Example_13_current_sense", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("13_current_sense_er95_n87.json");
    }
    TEST_CASE("Example_14_dab_xfmr", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("14_dab_xfmr_pm8770_n97.json");
    }
    TEST_CASE("Example_15_gan_inductor", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("15_gan_inductor_e138_3f4.json");
    }
    TEST_CASE("Example_16_coupled_inductor", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("16_coupled_inductor_e2513_dmr95.json");
    }
    TEST_CASE("Example_17_cllc_xfmr", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("17_cllc_xfmr_e5528_3c92a.json");
    }
    TEST_CASE("Example_18_stacked_inductor", "[example][magnetic-simulator][smoke-test]") {
        // ABT #1567: the example ran at 20 kHz, below N27's published loss data (25 kHz), which MKF
        // refuses (ABT #1456). The example now runs at 25 kHz, inside the fitted span, and must simulate.
        run_example_simulation("18_stacked_inductor_e7033_n27.json");
    }
    TEST_CASE("Example_19_multi_op_xfmr", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("19_multi_op_xfmr_etd3920_pc95.json");
    }
    TEST_CASE("Example_20_iso_buckboost", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("20_iso_buckboost_epc25_3c91.json");
    }
    TEST_CASE("Example_21_interleaved_flyback", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("21_interleaved_flyback_etd39_n87.json");
    }
    TEST_CASE("Example_22_margin_tape_forward", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("22_margin_tape_forward_e4218_3c95.json");
    }
    TEST_CASE("Example_23_interleaved_llc", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("23_interleaved_llc_pq3530_n97.json");
    }
    TEST_CASE("Example_24_margin_interleaved_flyback", "[example][magnetic-simulator][smoke-test]") {
        run_example_simulation("24_margin_interleaved_flyback_pq3230_3c94.json");
    }

    // ABT #838: the published core temperature counts the winding's heat, and every published
    // loss is evaluated at the temperature the COMBINED dissipation produces. Checked by
    // re-solving the same thermal network from the published losses, in both simulate modes:
    //   1. coreLosses.temperature is the hottest CORE node of the full network for those losses
    //      (not the whole part's hot spot, which is normally a turn, and not a core-only network);
    //   2. removing the winding's losses from that network cools the core: the published value
    //      is above the core-losses-only (P_cu = 0) core temperature;
    //   3. the core losses are the core-loss model evaluated at the published core temperature;
    //   4. the winding losses were evaluated at the network's mean turn temperature.
    TEST_CASE("Simulated core temperature includes the winding losses", "[processor][magnetic-simulator][smoke-test]") {
        const double tolerance = MagneticSimulator::kOperatingTemperatureTolerance;
        for (const std::string filename : {"01_simple_inductor_etd34_n87.json", "03_buck_inductor_pq3230_n95.json"}) {
            for (bool fastMode : {false, true}) {
                INFO(filename << " fastMode=" << fastMode);
                auto path = get_examples_dir() / filename;
                REQUIRE(std::filesystem::exists(path));
                auto mas = OpenMagneticsTesting::mas_loader(path.string());
                auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
                auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
                if (!magnetic.get_coil().get_turns_description()) {
                    magnetic.get_mutable_coil().wind();
                }

                MagneticSimulator simulator;
                auto simulated = simulator.simulate(inputs, magnetic, fastMode);
                REQUIRE(simulated.get_outputs().size() == inputs.get_operating_points().size());

                for (size_t opIdx = 0; opIdx < simulated.get_outputs().size(); ++opIdx) {
                    const auto& output = simulated.get_outputs()[opIdx];
                    const auto& operatingPoint = simulated.get_inputs().get_operating_points()[opIdx];
                    // The simulator must not hand back its internal core temperature as the ambient.
                    REQUIRE(operatingPoint.get_conditions().get_ambient_temperature() ==
                            inputs.get_operating_points()[opIdx].get_conditions().get_ambient_temperature());
                    const double publishedCoreTemperature = output.get_core_losses()->get_temperature();
                    const double windingLosses = output.get_winding_losses()->get_winding_losses();
                    REQUIRE(windingLosses > 0);

                    auto fullConfig = TemperatureConfig::fromSimulatedOutput(operatingPoint, output);
                    auto full = Temperature(magnetic, fullConfig).calculateTemperatures();
                    REQUIRE(full.maximumCoreTemperature);
                    REQUIRE(full.averageTurnTemperature);
                    CAPTURE(publishedCoreTemperature, full.maximumCoreTemperature.value(), full.maximumTemperature,
                            output.get_core_losses()->get_core_losses(), windingLosses);
                    // 1.
                    CHECK(std::fabs(publishedCoreTemperature - full.maximumCoreTemperature.value()) < tolerance);

                    // 2.
                    auto coreOnlyLossesConfig = fullConfig;
                    coreOnlyLossesConfig.windingLosses = 0;
                    coreOnlyLossesConfig.windingLossesOutput = std::nullopt;
                    auto withoutWinding = Temperature(magnetic, coreOnlyLossesConfig).calculateTemperatures();
                    REQUIRE(withoutWinding.maximumCoreTemperature);
                    CAPTURE(withoutWinding.maximumCoreTemperature.value());
                    CHECK(publishedCoreTemperature > withoutWinding.maximumCoreTemperature.value() + tolerance);

                    // 3.
                    auto excitation = operatingPoint.get_excitations_per_winding()[0];
                    REQUIRE(excitation.get_magnetic_flux_density());
                    double coreLossesAtPublishedTemperature = CoreLosses().calculate_core_losses(
                        magnetic.get_core(), excitation, publishedCoreTemperature).get_core_losses();
                    CHECK_THAT(output.get_core_losses()->get_core_losses(),
                               Catch::Matchers::WithinRel(coreLossesAtPublishedTemperature, 1e-9));

                    // 4.
                    REQUIRE(output.get_winding_losses()->get_temperature());
                    CAPTURE(output.get_winding_losses()->get_temperature().value(), full.averageTurnTemperature.value());
                    CHECK(std::fabs(output.get_winding_losses()->get_temperature().value() -
                                    full.averageTurnTemperature.value()) < tolerance);
                }
            }
        }
    }

    // ABT #838: the 2-argument calculate_core_losses (core advisers' loss filters, PyOM/WASM
    // calculate_core_losses) runs the SAME loss <-> temperature solve as simulate():
    //   - wound coil: identical core losses and core temperature to simulate(), above the
    //     core-losses-only core temperature (the winding's heat is counted);
    //   - coil with no turns: the core-only network, converged: re-solving that network from the
    //     returned losses gives back the returned temperature, and the losses are the model at it.
    TEST_CASE("calculate_core_losses solves the same loss-temperature loop as simulate", "[processor][magnetic-simulator][smoke-test]") {
        const double tolerance = MagneticSimulator::kOperatingTemperatureTolerance;
        // 05 (PFC inductor, T 40/20) replaces 03 (buck inductor, PQ 32/30): 03's core-only rise is
        // 0.6 K, inside the 2 K loop tolerance, so its core-only solve stops at ambient and the
        // "above ambient" check has nothing to measure. 05's core-only rise is ~15 K (Alf's decision).
        for (const std::string filename : {"01_simple_inductor_etd34_n87.json", "05_pfc_inductor_t4020_hf60.json"}) {
            INFO(filename);
            auto path = get_examples_dir() / filename;
            REQUIRE(std::filesystem::exists(path));
            auto mas = OpenMagneticsTesting::mas_loader(path.string());
            auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
            auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
            if (!magnetic.get_coil().get_turns_description()) {
                magnetic.get_mutable_coil().wind();
            }
            MagneticSimulator simulator;
            auto simulated = simulator.simulate(inputs, magnetic, true);
            const CoreLossesOutput simulatedCoreLosses = simulated.get_outputs()[0].get_core_losses().value();

            // Wound coil.
            auto operatingPoint = inputs.get_operating_points()[0];
            const double ambientTemperature = operatingPoint.get_conditions().get_ambient_temperature();
            auto wound = simulator.calculate_core_losses(operatingPoint, magnetic);
            CAPTURE(wound.get_temperature(), simulatedCoreLosses.get_temperature(), wound.get_core_losses(),
                    simulatedCoreLosses.get_core_losses());
            REQUIRE(operatingPoint.get_conditions().get_ambient_temperature() == ambientTemperature);
            CHECK_THAT(wound.get_temperature(), Catch::Matchers::WithinAbs(simulatedCoreLosses.get_temperature(), 1e-9));
            CHECK_THAT(wound.get_core_losses(), Catch::Matchers::WithinRel(simulatedCoreLosses.get_core_losses(), 1e-9));

            // Coil with no turns: core-only network, converged.
            OpenMagnetics::Magnetic unwoundMagnetic = magnetic;
            unwoundMagnetic.get_mutable_coil().set_turns_description(std::nullopt);
            unwoundMagnetic.get_mutable_coil().set_layers_description(std::nullopt);
            unwoundMagnetic.get_mutable_coil().set_sections_description(std::nullopt);
            auto unwoundOperatingPoint = inputs.get_operating_points()[0];
            auto coreOnly = simulator.calculate_core_losses(unwoundOperatingPoint, unwoundMagnetic);
            REQUIRE(unwoundOperatingPoint.get_conditions().get_ambient_temperature() == ambientTemperature);
            auto coreOnlyConfig = TemperatureConfig::fromMasOperatingConditions(unwoundOperatingPoint.get_conditions());
            coreOnlyConfig.plotSchematic = false;
            coreOnlyConfig.coreOnly = true;
            coreOnlyConfig.coreLosses = coreOnly.get_core_losses();
            auto coreOnlyNetwork = Temperature(unwoundMagnetic, coreOnlyConfig).calculateTemperatures();
            REQUIRE(coreOnlyNetwork.maximumCoreTemperature);
            CAPTURE(coreOnly.get_temperature(), coreOnlyNetwork.maximumCoreTemperature.value(), coreOnly.get_core_losses());
            CHECK(coreOnly.get_temperature() > ambientTemperature + tolerance);
            CHECK(std::fabs(coreOnly.get_temperature() - coreOnlyNetwork.maximumCoreTemperature.value()) < tolerance);
            auto excitation = unwoundOperatingPoint.get_excitations_per_winding()[0];
            REQUIRE(excitation.get_magnetic_flux_density());
            CHECK_THAT(coreOnly.get_core_losses(), Catch::Matchers::WithinRel(CoreLosses().calculate_core_losses(
                unwoundMagnetic.get_core(), excitation, coreOnly.get_temperature()).get_core_losses(), 1e-9));
            // The winding's heat warms the core.
            CHECK(wound.get_temperature() > coreOnly.get_temperature() + tolerance);
        }
    }

    TEST_CASE("The loss-temperature solve computes the winding's magnetic field once", "[processor][magnetic-simulator][smoke-test]") {
        // ABT #838: H at every turn follows from the currents and the geometry, not from the
        // temperature, so the solve computes it once and evaluates only the per-turn kernels at each
        // iterate's winding temperature. Counted through the field's turn-sum lookups: one solve must
        // make exactly the lookups of ONE winding-loss evaluation, although it iterates.
        const double tolerance = MagneticSimulator::kOperatingTemperatureTolerance;
        for (const std::string filename : {"01_simple_inductor_etd34_n87.json", "03_buck_inductor_pq3230_n95.json"}) {
            INFO(filename);
            auto path = get_examples_dir() / filename;
            REQUIRE(std::filesystem::exists(path));
            auto mas = OpenMagneticsTesting::mas_loader(path.string());
            auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
            auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
            if (!magnetic.get_coil().get_turns_description()) {
                magnetic.get_mutable_coil().wind();
            }
            MagneticSimulator simulator;
            auto lookups = []() {
                auto statistics = MagneticField::get_turn_sums_cache_statistics();
                return statistics.hits + statistics.misses;
            };

            auto operatingPoint = inputs.get_operating_points()[0];
            const double ambientTemperature = operatingPoint.get_conditions().get_ambient_temperature();
            size_t before = lookups();
            auto windingLossesAtAmbient = simulator.calculate_winding_losses(operatingPoint, magnetic, ambientTemperature);
            const size_t lookupsPerEvaluation = lookups() - before;
            REQUIRE(windingLossesAtAmbient.get_winding_losses() > 0);
            REQUIRE(lookupsPerEvaluation > 0);

            std::optional<WindingLossesOutput> solvedWindingLosses;
            before = lookups();
            auto solved = simulator.calculate_core_losses(operatingPoint, magnetic, &solvedWindingLosses);
            const size_t lookupsPerSolve = lookups() - before;
            REQUIRE(solvedWindingLosses);
            CAPTURE(solved.get_temperature(), lookupsPerEvaluation, lookupsPerSolve);
            // The first iterate is at ambient; a converged temperature this far from it took more
            // than one iterate, so a field computed per iterate would show up as extra lookups.
            REQUIRE(solved.get_temperature() > ambientTemperature + tolerance);
            CHECK(lookupsPerSolve == lookupsPerEvaluation);

            // The field AND the current split between a winding's parallels (ABT #1314) are those
            // of ambient at every iterate; only the conductors' resistance follows the iterate
            // (decided by Alf 2026-10-13). The solved winding temperature is read back from the
            // solved DC losses, which are linear in it.
            auto& settings = Settings::GetInstance();
            SettingsGuard<int> mirroringGuard(settings, &Settings::get_magnetic_field_mirroring_dimension,
                                              &Settings::set_magnetic_field_mirroring_dimension, 1);
            auto dcLosses = [&](double temperature) {
                double total = 0;
                const auto ohmicLosses = WindingOhmicLosses::calculate_ohmic_losses(magnetic.get_coil(), operatingPoint, temperature);
                const auto lossesPerWinding = ohmicLosses.get_winding_losses_per_winding().value();
                for (auto& winding : lossesPerWinding) {
                    total += winding.get_ohmic_losses()->get_losses();
                }
                return total;
            };
            double solvedDcLosses = 0;
            const auto solvedLossesPerWinding = solvedWindingLosses->get_winding_losses_per_winding().value();
            for (auto& winding : solvedLossesPerWinding) {
                solvedDcLosses += winding.get_ohmic_losses()->get_losses();
            }
            const double dcLossesAtAmbient = dcLosses(ambientTemperature);
            const double dcLossesPerKelvin = dcLosses(ambientTemperature + 100) / 100 - dcLossesAtAmbient / 100;
            const double solvedWindingTemperature = ambientTemperature + (solvedDcLosses - dcLossesAtAmbient) / dcLossesPerKelvin;
            REQUIRE(solvedWindingTemperature > ambientTemperature + tolerance);
            REQUIRE_THAT(dcLosses(solvedWindingTemperature), Catch::Matchers::WithinRel(solvedDcLosses, 1e-12));

            WindingLosses windingLossesModel;
            auto fieldAtAmbient = windingLossesModel.calculate_magnetic_field(magnetic, operatingPoint, ambientTemperature);
            const double lossesFieldAtAmbient = windingLossesModel.calculate_losses(magnetic, fieldAtAmbient, solvedWindingTemperature).get_winding_losses();
            const double lossesFieldAtSolved = windingLossesModel.calculate_losses(magnetic, operatingPoint, solvedWindingTemperature).get_winding_losses();
            CAPTURE(solvedWindingTemperature, solvedWindingLosses->get_winding_losses(), lossesFieldAtAmbient, lossesFieldAtSolved);
            CHECK_THAT(solvedWindingLosses->get_winding_losses(), Catch::Matchers::WithinRel(lossesFieldAtAmbient, 1e-9));
            bool hasParallels = false;
            for (auto& winding : magnetic.get_coil().get_functional_description()) {
                hasParallels = hasParallels || winding.get_number_parallels() > 1;
            }
            // 03_buck (12 turns x 2 parallels) must be able to tell the two choices apart.
            CHECK(hasParallels == (filename == "03_buck_inductor_pq3230_n95.json"));
            if (hasParallels) {
                CHECK(std::fabs(lossesFieldAtSolved - lossesFieldAtAmbient) > 1e-7 * lossesFieldAtAmbient);
                CHECK(std::fabs(lossesFieldAtSolved - lossesFieldAtAmbient) < 1e-4 * lossesFieldAtAmbient);
            }
        }
    }

    TEST_CASE("Inside a loss-solve memo scope an exact repeat returns the stored solve", "[processor][magnetic-simulator][smoke-test]") {
        // ABT #838: an adviser run solves the same candidate at the same operating point more than
        // once. Inside a LossSolveMemoScope the repeat returns the stored solve, which must be exactly
        // what a fresh solve returns; any difference in the inputs, the models or the Settings is a
        // different solve.
        auto path = get_examples_dir() / "01_simple_inductor_etd34_n87.json";
        REQUIRE(std::filesystem::exists(path));
        auto mas = OpenMagneticsTesting::mas_loader(path.string());
        auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
        auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
        if (!magnetic.get_coil().get_turns_description()) {
            magnetic.get_mutable_coil().wind();
        }
        const auto operatingPoint = inputs.get_operating_points()[0];
        const double ambientTemperature = operatingPoint.get_conditions().get_ambient_temperature();
        struct Solve {
            json coreLosses, windingLosses, temperature, excitation;
        };
        auto solve = [&magnetic](OperatingPoint thisOperatingPoint) {
            MagneticSimulator simulator;  // a fresh simulator each time, as different filters hold
            std::optional<WindingLossesOutput> windingLosses;
            std::optional<TemperatureOutput> temperature;
            auto coreLosses = simulator.calculate_core_losses(thisOperatingPoint, magnetic, &windingLosses, &temperature);
            REQUIRE(windingLosses);
            REQUIRE(temperature);
            Solve result;
            to_json(result.coreLosses, coreLosses);
            to_json(result.windingLosses, windingLosses.value());
            to_json(result.temperature, temperature.value());
            to_json(result.excitation, thisOperatingPoint.get_excitations_per_winding()[0]);
            return result;
        };
        auto same = [](const Solve& a, const Solve& b) {
            return a.coreLosses.dump() == b.coreLosses.dump() && a.windingLosses.dump() == b.windingLosses.dump() &&
                   a.temperature.dump() == b.temperature.dump() && a.excitation.dump() == b.excitation.dump();
        };

        REQUIRE(LossSolveMemoScope::get_statistics().entries == 0);
        auto fresh = solve(operatingPoint);
        {
            LossSolveMemoScope scope;
            auto first = solve(operatingPoint);
            auto repeat = solve(operatingPoint);
            auto statistics = LossSolveMemoScope::get_statistics();
            CHECK(statistics.misses == 1);
            CHECK(statistics.hits == 1);
            CHECK(same(first, fresh));
            CHECK(same(repeat, fresh));

            // A different operating point is a different solve.
            auto warmer = operatingPoint;
            warmer.get_mutable_conditions().set_ambient_temperature(ambientTemperature + 1);
            auto warmerSolve = solve(warmer);
            CHECK(!same(warmerSolve, fresh));
            statistics = LossSolveMemoScope::get_statistics();
            CHECK(statistics.misses == 2);
            CHECK(statistics.hits == 1);

            // So is a different setting.
            {
                auto& settings = Settings::GetInstance();
                SettingsGuard<bool> fringingGuard(settings, &Settings::get_magnetic_field_include_fringing,
                                                  &Settings::set_magnetic_field_include_fringing,
                                                  !settings.get_magnetic_field_include_fringing());
                solve(operatingPoint);
                statistics = LossSolveMemoScope::get_statistics();
                CHECK(statistics.misses == 3);
                CHECK(statistics.hits == 1);
            }

            // A nested scope uses the outer memo.
            {
                LossSolveMemoScope inner;
                auto nestedRepeat = solve(operatingPoint);
                CHECK(same(nestedRepeat, fresh));
                statistics = LossSolveMemoScope::get_statistics();
                CHECK(statistics.hits == 2);
                CHECK(statistics.entries == 3);
            }
        }
        // The memo ends with the outermost scope.
        CHECK(LossSolveMemoScope::get_statistics().entries == 0);
    }

    // =========================================================================
    // Datasheet builder: build_datasheet() turns a simulated MAS into a
    // manufacturer-style datasheetInfo block populated from the Outputs and the
    // magnetic geometry.
    // =========================================================================

    static OpenMagnetics::Mas simulate_example(const std::string& filename) {
        auto path = get_examples_dir() / filename;
        REQUIRE(std::filesystem::exists(path));

        auto mas = OpenMagneticsTesting::mas_loader(path.string());
        auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
        auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
        mas.set_magnetic(magnetic);
        mas.set_inputs(inputs);

        MagneticSimulator simulator;
        return simulator.simulate(mas.get_inputs(), mas.get_magnetic());
    }

    TEST_CASE("Build_Datasheet_Inductor", "[processor][magnetic-simulator][datasheet]") {
        auto mas = simulate_example("01_simple_inductor_etd34_n87.json");

        MagneticSimulator simulator;
        auto manufacturerInfo = simulator.build_datasheet(mas);

        // Attached back onto the magnetic.
        REQUIRE(mas.get_magnetic().get_manufacturer_info());
        REQUIRE(mas.get_magnetic().get_manufacturer_info()->get_datasheet_info());
        auto datasheet = mas.get_magnetic().get_manufacturer_info()->get_datasheet_info().value();

        // Part: single winding, material populated.
        REQUIRE(datasheet.get_part());
        REQUIRE(datasheet.get_part()->get_number_of_windings().value() == 1);
        REQUIRE(datasheet.get_part()->get_material());

        // Electrical: single-winding part uses dcResistance (not dcResistances), inductance > 0.
        REQUIRE(datasheet.get_electrical());
        auto electricalVec1 = datasheet.get_electrical().value();
        REQUIRE(!electricalVec1.empty());
        auto& electrical = electricalVec1[0];
        REQUIRE(electrical.get_inductance());
        REQUIRE(resolve_dimensional_values(electrical.get_inductance().value()) > 0);
        REQUIRE(electrical.get_dc_resistance());
        REQUIRE(resolve_dimensional_values(electrical.get_dc_resistance().value()) > 0);
        REQUIRE_FALSE(electrical.get_dc_resistances());
        REQUIRE(electrical.get_saturation_current_peak().value() > 0);

        // Rated (heat-limited) current and SRF / impedance sweep are populated.
        REQUIRE(electrical.get_rated_currents());
        REQUIRE(electrical.get_rated_currents()->at(0) > 0);
        REQUIRE(electrical.get_self_resonant_frequency().value() > 0);
        REQUIRE(electrical.get_impedance_points());
        REQUIRE(electrical.get_impedance_points()->size() > 1);
        REQUIRE(electrical.get_impedance_points()->front().get_impedance().get_magnitude() > 0);
        REQUIRE(electrical.get_maximum_impedance().value() > 0);

        // Mechanical: all three body dimensions in metres and positive.
        REQUIRE(datasheet.get_mechanical());
        auto mechanical = datasheet.get_mechanical().value();
        REQUIRE(resolve_dimensional_values(mechanical.get_width().value()) > 0);
        REQUIRE(resolve_dimensional_values(mechanical.get_height().value()) > 0);
        REQUIRE(resolve_dimensional_values(mechanical.get_length().value()) > 0);

        // Thermal: a temperature rise relative to ambient.
        REQUIRE(datasheet.get_thermal());
        REQUIRE(datasheet.get_thermal()->get_operating_temperature());
        REQUIRE(datasheet.get_thermal()->get_temperature_rise());

        // Application: switching frequency recorded.
        REQUIRE(datasheet.get_application());
        REQUIRE(datasheet.get_application()->get_switching_frequency().value() > 0);
    }

    TEST_CASE("Build_Datasheet_Transformer", "[processor][magnetic-simulator][datasheet]") {
        auto mas = simulate_example("02_flyback_efd25_3c95.json");

        MagneticSimulator simulator;
        simulator.build_datasheet(mas);

        auto datasheet = mas.get_magnetic().get_manufacturer_info()->get_datasheet_info().value();
        auto electricalVec = datasheet.get_electrical().value();
        REQUIRE(!electricalVec.empty());
        auto& electrical = electricalVec[0];

        // Multi-winding part: per-winding dcResistances, a turns ratio, no scalar dcResistance.
        REQUIRE(datasheet.get_part()->get_number_of_windings().value() > 1);
        REQUIRE(electrical.get_dc_resistances());
        REQUIRE(electrical.get_dc_resistances()->size() == static_cast<size_t>(datasheet.get_part()->get_number_of_windings().value()));
        REQUIRE_FALSE(electrical.get_dc_resistance());
        REQUIRE(electrical.get_turns_ratios());
    }

    TEST_CASE("Build_Datasheet_Throws_Without_Outputs", "[processor][magnetic-simulator][datasheet]") {
        auto path = get_examples_dir() / "01_simple_inductor_etd34_n87.json";
        auto mas = OpenMagneticsTesting::mas_loader(path.string());
        mas.set_outputs({});

        MagneticSimulator simulator;
        REQUIRE_THROWS(simulator.build_datasheet(mas));
    }

}

// ABT #366/#362/#357: the FULL simulate path over all four new families. Individual models were
// pinned family by family (inductance, saturation, losses, capacitance, leakage, impedance,
// temperature), but MagneticSimulator is what a real caller actually invokes — it runs the whole
// chain together and cross-feeds results (core losses -> temperature -> permeability -> losses).
// An integration that works model-by-model can still fall over here, and nothing had checked it.
TEST_CASE("Test_Magnetic_Simulator_New_Core_Families",
          "[processor][magnetic-simulator][drum][drum-ring][drum-semishielded][molded]") {
    settings.reset();
    clear_databases();

    auto buildCustomCore = [](json shapeJson, const std::string& coreType, json coating = json(),
                              const std::string& materialName = "3C90") {
        json coreJson;
        coreJson["functionalDescription"] = {
            {"type", coreType}, {"material", materialName}, {"shape", shapeJson},
            {"gapping", json::array()}, {"numberStacks", 1}};
        if (!coating.is_null()) {
            coreJson["functionalDescription"]["coating"] = coating;
        }
        OpenMagnetics::Core core(coreJson);
        core.process_data();
        core.process_gap();
        return core;
    };

    json drumDimensions = {
        {"A", {{"nominal", 0.0038}}}, {"B", {{"nominal", 0.0018}}}, {"C", {{"nominal", 0.0015}}},
        {"D", {{"nominal", 0.0004}}}, {"E", {{"nominal", 0.0010}}}, {"F", {{"nominal", 0.0004}}}};
    json semishieldedDimensions = drumDimensions;
    semishieldedDimensions["J"] = {{"nominal", 0.0040}};
    semishieldedDimensions["K"] = {{"nominal", 0.0040}};
    semishieldedDimensions["L"] = {{"nominal", 0.0018}};

    std::vector<std::pair<std::string, OpenMagnetics::Core>> cores;
    cores.emplace_back("drum", OpenMagneticsTesting::get_quick_core("DRH-14X20-4C", json::array(), 1, "3C90"));
    cores.emplace_back("drumRing", OpenMagneticsTesting::get_quick_core("DR 2.3 + SRI 3.0", json::array(), 1, "3C90"));
    cores.emplace_back("drumSemishielded", buildCustomCore(
        {{"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "drumSemishielded"},
         {"aliases", json::array()}, {"name", "LQS-like 4018"}, {"dimensions", semishieldedDimensions}},
        "pieceAndPlate", {{"type", "magneticEpoxy"}, {"thickness", {{"nominal", 0.0001}}}, {"material", "Kool M\u00b5 26"}}));
    // Molded bodies are a low-permeability metal composite, NOT ferrite: building this fixture
    // with a mu_i~2000 ferrite drove B (and hence core loss) an order of magnitude high for a
    // 4 mm part. Use a powder grade in the family's real permeability class instead.
    cores.emplace_back("molded", buildCustomCore(
        {{"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "molded"},
         {"aliases", json::array()}, {"name", "MAPI-like 4020"},
         {"dimensions", {
             {"A", {{"nominal", 0.0041}}}, {"B", {{"nominal", 0.0021}}}, {"C", {{"nominal", 0.0041}}},
             {"D", {{"nominal", 0.0014}}}, {"E", {{"nominal", 0.0030}}}, {"F", {{"nominal", 0.0012}}}}}},
        "closedShape", json(), "Kool M\u00b5 26"));

    for (auto& [label, core] : cores) {
        json coilJson;
        coilJson["bobbin"] = "Dummy";
        coilJson["functionalDescription"] = json::array({{
            {"name", "winding 0"}, {"numberTurns", 8}, {"numberParallels", 1},
            {"isolationSide", "primary"}, {"wire", "Round 0.1 - Grade 1"}}});
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(OpenMagnetics::Coil(coilJson, false));
        auto completed = OpenMagnetics::magnetic_autocomplete(magnetic);

        // ABT #1567: 400 kHz, inside 3C90's fitted Steinmetz span (25 kHz to 446.69 kHz); at the
        // former 500 kHz core losses are an extrapolation, which MKF refuses since ABT #1456.
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(
            400000, 5e-6, 25, WaveformLabel::TRIANGULAR, 0.2, 0.5, 0.2);
        auto completedInputs = OpenMagnetics::inputs_autocomplete(inputs, completed);

        OpenMagnetics::MagneticSimulator simulator;
        OpenMagnetics::Mas mas;
        REQUIRE_NOTHROW(mas = simulator.simulate(completedInputs, completed));
        REQUIRE(mas.get_outputs().size() > 0);
        auto& output = mas.get_outputs()[0];

        REQUIRE(output.get_inductance());
        double inductance = OpenMagnetics::resolve_dimensional_values(
            output.get_inductance()->get_magnetizing_inductance().get_magnetizing_inductance());
        REQUIRE(output.get_core_losses());
        double coreLosses = output.get_core_losses()->get_core_losses();
        REQUIRE(output.get_winding_losses());
        double windingLosses = output.get_winding_losses()->get_winding_losses();

        UNSCOPED_INFO(label << ": L " << inductance * 1e6 << " uH, core losses " << coreLosses * 1e3
                      << " mW, winding losses " << windingLosses * 1e3 << " mW");
        CHECK(std::isfinite(inductance));
        CHECK(inductance > 0);
        CHECK(std::isfinite(coreLosses));
        CHECK(coreLosses >= 0);
        CHECK(std::isfinite(windingLosses));
        CHECK(windingLosses > 0);
    }
    settings.reset();
}

// A web-builder design (buck, 98 P 11/9, 2x15 turns) whose core arrived WITHOUT a
// processedDescription — the builder cleared it while re-processing the core and a
// simulation fired in that window. A Magnetic deserialized from JSON does not process its
// core, and MagnetizingInductance read the effective area through the empty optional
// (UB): area 0, flux / 0, and "Waveform data contains NaN" far downstream. It must name
// the missing core data instead.
TEST_CASE("MagneticSimulator names an unprocessed core instead of dividing by a zero area",
          "[processor][magnetic-simulator][core-not-processed]") {
    std::ifstream fi(OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "web_buck_inputs_unprocessed_core_1369.json"));
    std::ifstream fm(OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "web_buck_magnetic_unprocessed_core_1369.json"));
    REQUIRE(fi.good());
    REQUIRE(fm.good());
    OpenMagnetics::Inputs inputs(nlohmann::json::parse(fi));
    OpenMagnetics::Magnetic magnetic(nlohmann::json::parse(fm));
    REQUIRE(!magnetic.get_core().get_processed_description());

    MagneticSimulator magneticSimulator;
    try {
        magneticSimulator.simulate(inputs, magnetic);
        FAIL("simulating a magnetic whose core has no processed description must throw");
    }
    catch (const CoreNotProcessedException& e) {
        CHECK(std::string(e.what()).find("processedDescription") != std::string::npos);
    }
}

// ABT #1855 (Alf, 2026-10-14): a saturating operating point gets its numbers and a saturation finding
// for THAT point; it no longer stops the simulation of the others. The MAS example flyback (EFD 25, 3C95),
// whose operating point simulates cool (about 69 C, B_pk 0.28 T), with a second operating point that is
// the first with twice its winding currents (the flyback's magnetizing current is its winding currents'):
// simulate used to throw for the whole design.
TEST_CASE("Simulate gives every operating point its numbers and analyze_results names the saturating one",
          "[processor][magnetic-simulator][abt-1855]") {
    settings.reset();
    auto path = std::filesystem::path{std::source_location::current().file_name()}
                    .parent_path().append("..").append("MAS").append("examples").append("02_flyback_efd25_3c95.json");
    REQUIRE(std::filesystem::exists(path));
    auto mas = OpenMagneticsTesting::mas_loader(path.string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    json inputsJson;
    to_json(inputsJson, mas.get_inputs());
    REQUIRE(inputsJson["operatingPoints"].size() == 1);
    json saturatingJson = inputsJson["operatingPoints"][0];
    saturatingJson["name"] = "Twice the winding currents";
    for (auto& excitation : saturatingJson["excitationsPerWinding"]) {
        for (auto& value : excitation["current"]["waveform"]["data"]) {
            value = 2.0 * value.get<double>();
        }
        excitation["current"].erase("processed");
        excitation["current"].erase("harmonics");
        excitation.erase("magnetizingCurrent");
        excitation.erase("magneticFluxDensity");
    }
    inputsJson["operatingPoints"].push_back(saturatingJson);
    auto inputs = OpenMagnetics::inputs_autocomplete(OpenMagnetics::Inputs(inputsJson), magnetic);
    REQUIRE(inputs.get_operating_points().size() == 2);

    MagneticSimulator magneticSimulator;
    OpenMagnetics::Mas simulated;
    REQUIRE_NOTHROW(simulated = magneticSimulator.simulate(inputs, magnetic));
    REQUIRE(simulated.get_outputs().size() == 2);
    for (size_t index = 0; index < 2; ++index) {
        INFO("operating point " << index);
        auto solvedExcitation = simulated.get_inputs().get_operating_points()[index].get_excitations_per_winding()[0];
        UNSCOPED_INFO("operating point " << index << ": core losses at " << simulated.get_outputs()[index].get_core_losses()->get_temperature()
             << " C, Bpk " << (solvedExcitation.get_magnetic_flux_density() ? solvedExcitation.get_magnetic_flux_density()->get_processed()->get_peak().value() : -1.0) << " T");
        REQUIRE(simulated.get_outputs()[index].get_core_losses());
        CHECK(std::isfinite(simulated.get_outputs()[index].get_core_losses()->get_core_losses()));
        REQUIRE(simulated.get_outputs()[index].get_temperature());
    }

    auto findings = ResultsAnalyzer::analyze_results(simulated);
    std::vector<ResultFinding> saturationFindings;
    std::string allFindings;
    for (const auto& finding : findings) {
        allFindings += "\n" + to_string(finding.code) + ": " + finding.message;
        if (finding.code == ResultFindingCode::CORE_SATURATES) {
            saturationFindings.push_back(finding);
        }
    }
    INFO("findings:" << allFindings);
    REQUIRE(saturationFindings.size() == 1);
    CHECK(saturationFindings[0].operatingPointIndex == 1);
    REQUIRE(saturationFindings[0].operatingPointName);
    CHECK(saturationFindings[0].operatingPointName.value() == "Twice the winding currents");
    CHECK_THAT(saturationFindings[0].message, Catch::Matchers::ContainsSubstring("The core saturates") &&
                                              Catch::Matchers::ContainsSubstring("operating point 1"));

    // As JSON, the web's form: one object per finding naming its operating point.
    json findingsJson;
    to_json(findingsJson, findings);
    REQUIRE(findingsJson.is_array());
    REQUIRE(findingsJson.size() == findings.size());
    for (const auto& entry : findingsJson) {
        CHECK(entry.contains("operatingPointIndex"));
        CHECK(entry.contains("code"));
        CHECK(entry.contains("message"));
    }
    settings.reset();
}

// ABT #2109 (Alf, 2026-10-27): since the overlap reporting of ABT #2071, analyze_results threw for a
// stored design whose turns store no cross-sectional shape ("Turn 'Primary parallel 0 turn 0' has no
// cross-sectional shape"), and every finding was lost. A turn without a stored cross-section takes its
// outline from its winding's wire; the stored simple_inductor (round wire, outer diameter 0.676 mm,
// turns without crossSectionalShape) returns its findings again.
TEST_CASE("analyze_results judges a stored design whose turns store no cross-section by its wire's outline",
          "[processor][results-analyzer][abt-2109]") {
    settings.reset();
    auto path = std::filesystem::path{std::source_location::current().file_name()}
                    .parent_path().append("testData").append("simple_inductor.json");
    REQUIRE(std::filesystem::exists(path));
    auto mas = OpenMagneticsTesting::mas_loader(path.string());
    const auto storedTurns = mas.get_magnetic().get_coil().get_turns_description().value();
    REQUIRE(storedTurns.size() == 31);
    for (const auto& turn : storedTurns) {
        REQUIRE_FALSE(turn.get_cross_sectional_shape());   // the premise: the stored turns carry no shape
    }
    auto coil = mas.get_magnetic().get_coil();
    CHECK_FALSE(coil.turn_outlines_unknown_reason());
    // As the web does: simulate the stored design, then analyze it.
    MagneticSimulator magneticSimulator;
    OpenMagnetics::Mas simulated;
    REQUIRE_NOTHROW(simulated = magneticSimulator.simulate(mas.get_inputs(), mas.get_magnetic()));
    const auto simulatedTurns = simulated.get_magnetic().get_coil().get_turns_description().value();   // bound
    REQUIRE(simulatedTurns.size() == 31);
    for (const auto& turn : simulatedTurns) {
        REQUIRE_FALSE(turn.get_cross_sectional_shape());   // simulate keeps the stored turns
    }
    std::vector<ResultFinding> findings;
    REQUIRE_NOTHROW(findings = ResultsAnalyzer::analyze_results(simulated));
    std::string allFindings;
    std::set<ResultFindingCode> codes;
    for (const auto& finding : findings) {
        allFindings += "\n" + to_string(finding.code) + ": " + finding.message;
        codes.insert(finding.code);
    }
    INFO("findings:" << allFindings);
    CHECK(codes.count(ResultFindingCode::CORE_SATURATES) == 1);
    CHECK(codes.count(ResultFindingCode::CORE_ABOVE_MATERIAL_DATA_LIMIT) == 1);
    // Its wire states no temperature limit (no coating rating, no rated material), so the 543 C winding
    // is reported as unchecked rather than above a limit.
    CHECK(codes.count(ResultFindingCode::WINDING_TEMPERATURE_LIMIT_NOT_STATED) == 1);
    CHECK(codes.count(ResultFindingCode::TURN_OUTLINES_NOT_EVALUABLE) == 0);
    settings.reset();
}

// ABT #2109 (Alf, 2026-10-27): when the wire does not give the outline either, the overlap and
// outside-window checks are not evaluable for the design: ONE finding says so, with the reason, and
// every other finding is still returned. simple_inductor with its round wire's outer diameter removed.
TEST_CASE("analyze_results reports the turn overlaps not evaluable when neither the turns nor their wire give an outline, and keeps every other finding",
          "[processor][results-analyzer][abt-2109]") {
    settings.reset();
    auto path = std::filesystem::path{std::source_location::current().file_name()}
                    .parent_path().append("testData").append("simple_inductor.json");
    REQUIRE(std::filesystem::exists(path));
    auto mas = OpenMagneticsTesting::mas_loader(path.string());
    auto magnetic = mas.get_magnetic();
    auto coil = magnetic.get_coil();
    auto wire = coil.resolve_wire(size_t(0));
    REQUIRE(wire.get_type() == WireType::ROUND);
    REQUIRE(wire.get_outer_diameter());
    wire.set_outer_diameter(std::nullopt);
    coil.get_mutable_functional_description()[0].set_wire(wire);
    magnetic.set_coil(coil);
    mas.set_magnetic(magnetic);
    REQUIRE_FALSE(mas.get_magnetic().get_coil().get_turns_description().value()[0].get_cross_sectional_shape());

    auto judged = mas.get_magnetic().get_coil();
    const auto reason = judged.turn_outlines_unknown_reason();
    REQUIRE(reason);
    CHECK_THAT(reason.value(), Catch::Matchers::ContainsSubstring("Primary parallel 0 turn 0") &&
                               Catch::Matchers::ContainsSubstring("no outer diameter"));
    // The judges themselves still throw: only analyze_results turns it into a finding.
    CHECK_THROWS_AS(judged.judge_copper_overlaps(), InvalidInputException);

    MagneticSimulator magneticSimulator;
    OpenMagnetics::Mas simulated;
    REQUIRE_NOTHROW(simulated = magneticSimulator.simulate(mas.get_inputs(), mas.get_magnetic()));
    REQUIRE_FALSE(simulated.get_magnetic().get_coil().get_turns_description().value()[0].get_cross_sectional_shape());
    std::vector<ResultFinding> findings;
    REQUIRE_NOTHROW(findings = ResultsAnalyzer::analyze_results(simulated));
    std::string allFindings;
    std::vector<ResultFinding> notEvaluable;
    std::set<ResultFindingCode> codes;
    for (const auto& finding : findings) {
        allFindings += "\n" + to_string(finding.code) + ": " + finding.message;
        codes.insert(finding.code);
        if (finding.code == ResultFindingCode::TURN_OUTLINES_NOT_EVALUABLE) {
            notEvaluable.push_back(finding);
        }
    }
    INFO("findings:" << allFindings);
    REQUIRE(notEvaluable.size() == 1);
    CHECK_FALSE(notEvaluable[0].operatingPointIndex);
    CHECK_THAT(notEvaluable[0].message, Catch::Matchers::ContainsSubstring("cannot be checked") &&
                                        Catch::Matchers::ContainsSubstring("Primary parallel 0 turn 0") &&
                                        Catch::Matchers::ContainsSubstring("no outer diameter") &&
                                        !Catch::Matchers::ContainsSubstring("ABT"));
    CHECK(codes.count(ResultFindingCode::CORE_SATURATES) == 1);
    CHECK(codes.count(ResultFindingCode::CORE_ABOVE_MATERIAL_DATA_LIMIT) == 1);
    // Its wire states no temperature limit (no coating rating, no rated material), so the 543 C winding
    // is reported as unchecked rather than above a limit.
    CHECK(codes.count(ResultFindingCode::WINDING_TEMPERATURE_LIMIT_NOT_STATED) == 1);
    CHECK(codes.count(ResultFindingCode::TURNS_OVERLAP) == 0);
    CHECK(codes.count(ResultFindingCode::TURNS_OUTSIDE_WINDOW) == 0);
    settings.reset();
}

// ABT #1855, second saturation criterion (Alf, 2026-10-14): a DC bias that drives the material's
// permeability model to or below the vacuum's (main's SaturationException condition until #1855) is a
// saturation even when the flux density on the material's magnetisation curve stays under B_sat. 3C95
// on an ER 28 with a 1 mm gap and 12 turns: the model takes mu_r = 1 from ~44.6 A, while B stays under
// 0.53 T up to 260 A, so the B criterion alone would let the air-core point pass.
TEST_CASE("A DC bias that collapses the permeability to mu_r = 1 is a saturation with B under B_sat",
          "[processor][magnetic-simulator][abt-1855]") {
    settings.reset();
    std::string shapeName = "ER 28";
    auto coil = OpenMagneticsTesting::get_quick_coil({12}, {1}, shapeName, 1, WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING, CoilAlignment::SPREAD,
                                                     CoilAlignment::CENTERED, {find_wire_by_name("Round 0.25 - FIW 6")}, true);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, OpenMagneticsTesting::get_ground_gap(0.001), 1, "3C95");
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    auto coreMaterial = core.resolve_material();
    const double temperature = 25;
    const double frequency = 100000;

    MagnetizingInductance model;
    double staticInductance = resolve_dimensional_values(model.calculate_inductance_from_number_turns_and_gapping(core, coil).get_magnetizing_inductance());
    auto evaluate = [&](double currentOffset) {
        auto operatingPoint = OpenMagnetics::Inputs::create_operating_point_with_sinusoidal_current_mask(
            frequency, staticInductance, temperature, {}, {sqrt(2)}, currentOffset);
        auto fluxDensity = model.calculate_inductance_and_magnetic_flux_density(core, coil, &operatingPoint).second;
        return std::make_tuple(operatingPoint, fluxDensity, model.get_last_dc_bias_point());
    };

    SECTION("check_saturation: B alone passes the collapsed point, the DC bias point flags it") {
        auto [operatingPoint, fluxDensity, dcBiasPoint] = evaluate(100);
        REQUIRE(dcBiasPoint);
        double peak = ResultsAnalyzer::get_peak_flux_density(fluxDensity);
        INFO("B peak " << peak << " T at 100 A, H_dc " << dcBiasPoint->magneticFieldDcBias << " A/m");
        auto byFluxDensity = ResultsAnalyzer::check_saturation(coreMaterial, peak, temperature);
        CHECK_FALSE(byFluxDensity.saturates);
        auto verdict = ResultsAnalyzer::check_saturation(coreMaterial, peak, temperature, 1.0, dcBiasPoint);
        INFO(verdict.message);
        CHECK(verdict.saturates);
        CHECK(verdict.permeabilityCollapsed);
        CHECK_FALSE(verdict.fluxDensityAboveSaturation);
        REQUIRE(verdict.biasedPermeability);
        CHECK(verdict.biasedPermeability.value() <= 1);
        CHECK_THAT(verdict.message, Catch::Matchers::ContainsSubstring("its permeability collapses"));
        // The model's own answer there is the air core.
        double airCoreInductance = 144 / ReluctanceModel::factory()->get_core_reluctance(core, 1.0).get_core_reluctance();
        auto output = model.calculate_inductance_from_number_turns_and_gapping(core, coil, &operatingPoint);
        CHECK_THAT(resolve_dimensional_values(output.get_magnetizing_inductance()), Catch::Matchers::WithinRel(airCoreInductance, 1e-9));
    }

    SECTION("analyze_saturation reads the DC point itself and names the operating point; an unsaturated bias has no finding") {
        for (double currentOffset : {100.0, 5.0}) {
            auto [operatingPoint, fluxDensity, dcBiasPoint] = evaluate(currentOffset);
            CoreLossesOutput coreLosses;
            coreLosses.set_core_losses(0);
            coreLosses.set_method_used("test");
            coreLosses.set_origin(ResultOrigin::SIMULATION);
            coreLosses.set_temperature(temperature);
            coreLosses.set_magnetic_flux_density(fluxDensity);
            OpenMagnetics::Outputs output;
            output.set_core_losses(coreLosses);
            auto findings = ResultsAnalyzer::analyze_saturation(magnetic, operatingPoint, output, 3);
            INFO("offset " << currentOffset << " A");
            if (currentOffset > 50) {
                REQUIRE(findings.size() == 1);
                CHECK(findings[0].code == ResultFindingCode::CORE_SATURATES);
                CHECK(findings[0].operatingPointIndex == 3);
                CHECK_THAT(findings[0].message, Catch::Matchers::ContainsSubstring("its permeability collapses"));
            }
            else {
                CHECK(findings.empty());
            }
        }
    }
    settings.reset();
}

// ABT #1876 (Alf decision 2026-10-14): an operating point with no thermal operating point (the thermal
// network does not converge on its losses inside the loss loop, or its winding runs away) affects only
// that point: its core and winding losses are those at ambient (as before ABT #838), it gets no
// temperature output, and analyze_results reports NO_THERMAL_OPERATING_POINT for it with MKF's message.
// The other points simulate normally; any other error still refuses the whole design.
namespace {
    // The WebFrontend's default design (tests/fixtures/default_two_winding_pq2715.json there): PQ 27/15,
    // 13 + 13 turns of the "Dummy" wire, +-5 A triangular at 100 kHz, 25 C. About 2 kW of copper.
    std::pair<OpenMagnetics::Magnetic, OpenMagnetics::Inputs> load_web_default_design() {
        auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "default_two_winding_pq2715.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);
        auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
        auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
        return {magnetic, inputs};
    }

    json read_web_default_fixture() {
        auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "default_two_winding_pq2715.json");
        std::ifstream file(path);
        return json::parse(file);
    }

    // The fixture's operating point as stored, with every winding current and voltage scaled and the derived signals
    // dropped: the Inputs constructor and inputs_autocomplete derive them again.
    json scaled_operating_point(double factor, const std::string& name) {
        json operatingPointJson = read_web_default_fixture()["inputs"]["operatingPoints"][0];
        operatingPointJson["name"] = name;
        for (auto& excitation : operatingPointJson["excitationsPerWinding"]) {
            for (const std::string signal : {"current", "voltage"}) {
                if (!excitation.contains(signal)) {
                    continue;
                }
                for (auto& value : excitation[signal]["waveform"]["data"]) {
                    value = factor * value.get<double>();
                }
                excitation[signal].erase("processed");
                excitation[signal].erase("harmonics");
            }
            excitation.erase("magnetizingCurrent");
            excitation.erase("magneticFluxDensity");
            excitation.erase("magneticFieldStrength");
        }
        return operatingPointJson;
    }

    // The fixture's inputs with these operating points, processed as mas_loader processes the fixture's own
    // (the fixture has no outputs: Inputs with the magnetizing inductance of its turns and gapping), then autocompleted.
    OpenMagnetics::Inputs web_default_inputs_with(const json& operatingPoints, const OpenMagnetics::Magnetic& magnetic) {
        json inputsJson = read_web_default_fixture()["inputs"];
        inputsJson["operatingPoints"] = operatingPoints;
        double magnetizingInductance = resolve_dimensional_values(MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(
            magnetic.get_core(), magnetic.get_coil()).get_magnetizing_inductance());
        return OpenMagnetics::inputs_autocomplete(OpenMagnetics::Inputs(inputsJson, true, magnetizingInductance), magnetic);
    }

    std::string describe_findings(const std::vector<ResultFinding>& findings) {
        std::string all;
        for (const auto& finding : findings) {
            all += "\n[" + (finding.operatingPointIndex ? std::to_string(finding.operatingPointIndex.value()) : std::string("magnetic")) + "] " + to_string(finding.code) + ": " + finding.message;
        }
        return all;
    }
}

TEST_CASE("The web's default design simulates with its losses at ambient and a NO_THERMAL_OPERATING_POINT finding",
          "[processor][magnetic-simulator][abt-1876]") {
    // Alf decision 2026-10-14, ABT #1876. On 98e62c64 simulate threw for the whole design: "Thermal network
    // did not converge in 100 iterations ... Loss loop at iterate 2 (core losses at 122.5 C, winding losses at
    // 13518.7 C)". Prod 426ec505 (losses at ambient) gave 1996.6 W of copper.
    settings.reset();
    auto [magnetic, inputs] = load_web_default_design();
    REQUIRE(inputs.get_operating_points().size() == 1);
    const double ambientTemperature = inputs.get_operating_points()[0].get_conditions().get_ambient_temperature();
    REQUIRE(ambientTemperature == 25);

    MagneticSimulator magneticSimulator;
    OpenMagnetics::Mas simulated;
    REQUIRE_NOTHROW(simulated = magneticSimulator.simulate(inputs, magnetic));
    REQUIRE(simulated.get_outputs().size() == 1);
    const auto& output = simulated.get_outputs()[0];
    CHECK_FALSE(output.get_temperature());
    REQUIRE(output.get_core_losses());
    REQUIRE(output.get_winding_losses());
    CHECK(output.get_core_losses()->get_temperature() == ambientTemperature);

    // The same losses evaluated at ambient, independently of the loss loop.
    OpenMagnetics::Magnetic woundMagnetic = magnetic;
    if (!woundMagnetic.get_coil().get_turns_description()) {
        woundMagnetic.get_mutable_coil().wind();
    }
    OperatingPoint operatingPoint = inputs.get_operating_points()[0];
    double windingLossesAtAmbient = MagneticSimulator().calculate_winding_losses(operatingPoint, woundMagnetic, ambientTemperature).get_winding_losses();
    MagnetizingInductance magnetizingInductanceModel(to_string(settings.get_reluctance_model()));
    magnetizingInductanceModel.set_topology(inputs.get_design_requirements().get_topology());
    OperatingPoint operatingPointForCore = inputs.get_operating_points()[0];
    auto fluxDensity = magnetizingInductanceModel.calculate_inductance_and_magnetic_flux_density(woundMagnetic, &operatingPointForCore).second;
    auto excitation = operatingPointForCore.get_excitations_per_winding()[0];
    excitation.set_magnetic_flux_density(fluxDensity);
    double coreLossesAtAmbient = CoreLosses().calculate_core_losses(woundMagnetic.get_core(), excitation, ambientTemperature).get_core_losses();
    INFO("winding " << output.get_winding_losses()->get_winding_losses() << " W (at ambient " << windingLossesAtAmbient
         << " W), core " << output.get_core_losses()->get_core_losses() << " W (at ambient " << coreLossesAtAmbient << " W)");
    CHECK_THAT(output.get_winding_losses()->get_winding_losses(), Catch::Matchers::WithinRel(windingLossesAtAmbient, 1e-12));
    CHECK_THAT(output.get_core_losses()->get_core_losses(), Catch::Matchers::WithinRel(coreLossesAtAmbient, 1e-12));
    // Prod 426ec505: 1996.6 W of copper.
    CHECK_THAT(output.get_winding_losses()->get_winding_losses(), Catch::Matchers::WithinRel(1996.6, 0.01));

    auto findings = ResultsAnalyzer::analyze_results(simulated);
    INFO("findings:" << describe_findings(findings));
    std::vector<ResultFinding> noThermal;
    for (const auto& finding : findings) {
        if (finding.code == ResultFindingCode::NO_THERMAL_OPERATING_POINT) {
            noThermal.push_back(finding);
        }
        // No temperature output: no temperature verdicts for it.
        CHECK(finding.code != ResultFindingCode::WINDING_TEMPERATURE_ABOVE_LIMIT);
        CHECK(finding.code != ResultFindingCode::WINDING_TEMPERATURE_LIMIT_NOT_STATED);
    }
    REQUIRE(noThermal.size() == 1);
    CHECK(noThermal[0].operatingPointIndex == 0);
    CHECK_THAT(noThermal[0].message, Catch::Matchers::ContainsSubstring("No steady temperature") &&
                                     Catch::Matchers::ContainsSubstring("the thermal model does not settle") &&
                                     Catch::Matchers::ContainsSubstring("at the ambient temperature"));
    REQUIRE(noThermal[0].details);
    CHECK_THAT(noThermal[0].details.value(), Catch::Matchers::ContainsSubstring("History:") &&
                                             Catch::Matchers::ContainsSubstring("Pcu="));
    json findingsJson;
    to_json(findingsJson, findings);
    bool codeInJson = false;
    for (const auto& entry : findingsJson) {
        codeInJson = codeInJson || entry["code"] == "NO_THERMAL_OPERATING_POINT";
    }
    CHECK(codeInJson);
    settings.reset();
}

TEST_CASE("analyze_results reads why a point has no thermal operating point from simulate's output, without solving it again",
          "[processor][magnetic-simulator][abt-2075]") {
    // ABT #2075 (Alf, 2026-10-26): simulate records the reason on the point's core losses (methodUsed: the loss model's,
    // then MagneticSimulator::kLossesAtAmbientMethodNote, then MKF's message), and analyze_results quotes it. A second
    // solve doubled the cost of every such point and could disagree with simulate (ABT #2060).
    settings.reset();
    auto [magnetic, inputs] = load_web_default_design();
    OpenMagnetics::Mas simulated;
    REQUIRE_NOTHROW(simulated = MagneticSimulator().simulate(inputs, magnetic));
    REQUIRE(simulated.get_outputs().size() == 1);
    REQUIRE_FALSE(simulated.get_outputs()[0].get_temperature());
    REQUIRE(simulated.get_outputs()[0].get_core_losses());
    auto coreLosses = simulated.get_outputs()[0].get_core_losses().value();
    INFO("methodUsed: " << coreLosses.get_method_used());
    auto reason = MagneticSimulator::get_no_thermal_operating_point_reason(coreLosses);
    REQUIRE(reason);
    CHECK_THAT(reason.value(), Catch::Matchers::StartsWith("No steady temperature") &&
                               Catch::Matchers::ContainsSubstring("the thermal model does not settle"));
    // The loss model's own name comes first.
    const size_t note = coreLosses.get_method_used().find(MagneticSimulator::kLossesAtAmbientMethodNote);
    REQUIRE(note != std::string::npos);
    CHECK(note > 0);

    // The finding quotes what the output records: a solve of its own would carry the solver's message instead.
    const std::string recorded = "No thermal operating point: the reason this output records.";
    coreLosses.set_method_used(coreLosses.get_method_used().substr(0, note) + MagneticSimulator::kLossesAtAmbientMethodNote + recorded);
    auto outputs = simulated.get_outputs();
    outputs[0].set_core_losses(coreLosses);
    simulated.set_outputs(outputs);
    auto findings = ResultsAnalyzer::analyze_results(simulated);
    INFO("findings:" << describe_findings(findings));
    size_t noThermal = 0;
    for (const auto& finding : findings) {
        if (finding.code == ResultFindingCode::NO_THERMAL_OPERATING_POINT) {
            ++noThermal;
            CHECK_THAT(finding.message, Catch::Matchers::ContainsSubstring(recorded) &&
                                        !Catch::Matchers::ContainsSubstring("the thermal model does not settle"));
        }
    }
    CHECK(noThermal == 1);
    settings.reset();
}

TEST_CASE("A point with no thermal operating point leaves the other points as they simulate alone",
          "[processor][magnetic-simulator][abt-1876]") {
    // Alf decision 2026-10-14, ABT #1876: point 0 is the web default with a hundredth of its currents (0.23 W of copper, it has a
    // thermal operating point at about 58 C), point 1 the web default itself (it has none).
    settings.reset();
    auto magnetic = load_web_default_design().first;
    json coolPoint = scaled_operating_point(0.01, "A hundredth of the current");
    json hotPoint = scaled_operating_point(1.0, "The web default");
    auto twoPoints = web_default_inputs_with(json::array({coolPoint, hotPoint}), magnetic);
    auto onePoint = web_default_inputs_with(json::array({coolPoint}), magnetic);

    SECTION("point 0 is bit-identical to simulating it alone; point 1 has ambient losses, no temperature and the finding") {
        OpenMagnetics::Mas alone;
        OpenMagnetics::Mas both;
        REQUIRE_NOTHROW(alone = MagneticSimulator().simulate(onePoint, magnetic));
        REQUIRE_NOTHROW(both = MagneticSimulator().simulate(twoPoints, magnetic));
        REQUIRE(alone.get_outputs().size() == 1);
        REQUIRE(both.get_outputs().size() == 2);
        REQUIRE(alone.get_outputs()[0].get_temperature());
        INFO("point 0: core " << alone.get_outputs()[0].get_core_losses()->get_core_losses() << " W, winding "
             << alone.get_outputs()[0].get_winding_losses()->get_winding_losses() << " W, hottest "
             << alone.get_outputs()[0].get_temperature()->get_maximum_temperature() << " C");
        json aloneOutput, bothOutput, aloneOperatingPoint, bothOperatingPoint;
        to_json(aloneOutput, alone.get_outputs()[0]);
        to_json(bothOutput, both.get_outputs()[0]);
        to_json(aloneOperatingPoint, alone.get_inputs().get_operating_points()[0]);
        to_json(bothOperatingPoint, both.get_inputs().get_operating_points()[0]);
        CHECK(aloneOutput.dump() == bothOutput.dump());
        CHECK(aloneOperatingPoint.dump() == bothOperatingPoint.dump());

        CHECK_FALSE(both.get_outputs()[1].get_temperature());
        CHECK(both.get_outputs()[1].get_core_losses()->get_temperature() ==
              both.get_inputs().get_operating_points()[1].get_conditions().get_ambient_temperature());

        auto findings = ResultsAnalyzer::analyze_results(both);
        INFO("findings:" << describe_findings(findings));
        size_t noThermalAt1 = 0;
        for (const auto& finding : findings) {
            CHECK_FALSE((finding.code == ResultFindingCode::NO_THERMAL_OPERATING_POINT && finding.operatingPointIndex == 0));
            if (finding.code == ResultFindingCode::NO_THERMAL_OPERATING_POINT && finding.operatingPointIndex == 1) {
                ++noThermalAt1;
                REQUIRE(finding.operatingPointName);
                CHECK(finding.operatingPointName.value() == "The web default");
            }
        }
        CHECK(noThermalAt1 == 1);
    }

    SECTION("fast mode publishes no temperature either, and analyze_results does not take it for no thermal operating point") {
        OpenMagnetics::Mas fast;
        REQUIRE_NOTHROW(fast = MagneticSimulator().simulate(twoPoints, magnetic, true));
        REQUIRE(fast.get_outputs().size() == 2);
        CHECK_FALSE(fast.get_outputs()[0].get_temperature());
        CHECK_FALSE(fast.get_outputs()[1].get_temperature());
        try {
            ResultsAnalyzer::analyze_results(fast);
            FAIL("analyze_results read a fast-mode output as one with no thermal operating point");
        }
        catch (const OpenMagnetics::InvalidInputException& error) {
            CHECK_THAT(std::string(error.what()), Catch::Matchers::ContainsSubstring("fast mode") &&
                                                  Catch::Matchers::ContainsSubstring("operating point 'A hundredth of the current'"));
        }
    }
    settings.reset();
}

TEST_CASE("A malformed operating point still refuses the whole design next to one with no thermal operating point",
          "[processor][magnetic-simulator][abt-1876]") {
    // Alf decision 2026-10-14, ABT #1876: only thermal non-convergence and the winding runaway are per point.
    // Point 0 has no thermal operating point (the web default); point 1 has no current in its second winding.
    // simulate must throw point 1's error, not stop at point 0's thermal one.
    settings.reset();
    auto magnetic = load_web_default_design().first;
    json hotPoint = scaled_operating_point(1.0, "The web default");
    auto twoPoints = web_default_inputs_with(json::array({hotPoint, hotPoint}), magnetic);
    twoPoints.get_mutable_operating_points()[1].get_mutable_excitations_per_winding()[1].set_current(std::nullopt);

    try {
        MagneticSimulator().simulate(twoPoints, magnetic);
        FAIL("simulate accepted a winding without a current");
    }
    catch (const OpenMagnetics::NoThermalOperatingPointException& error) {
        FAIL("simulate stopped at point 0's thermal error instead of point 1's malformed input: " << error.what());
    }
    catch (const OpenMagnetics::OpenMagneticsException& error) {
        INFO(error.what());
        CHECK_THAT(std::string(error.what()), !Catch::Matchers::ContainsSubstring("the thermal model does not settle"));
        CHECK_THAT(std::string(error.what()), Catch::Matchers::ContainsSubstring("urrent"));
    }
    settings.reset();
}

TEST_CASE("calculate_core_losses of a point with no thermal operating point returns its core losses at ambient and no temperature",
          "[processor][magnetic-simulator][abt-1876]") {
    // Alf decision 2026-10-14, ABT #1876: the public per-point entry point answers as simulate does per point. On 28d80298
    // the standalone call (the WASM and PyOM calculate_core_losses bindings) still threw "[CALCULATION_INVALID_RESULT]
    // No thermal operating point: ..." for the web's default design.
    settings.reset();
    auto [magnetic, inputs] = load_web_default_design();
    const double ambientTemperature = inputs.get_operating_points()[0].get_conditions().get_ambient_temperature();
    REQUIRE(ambientTemperature == 25);
    OpenMagnetics::Magnetic woundMagnetic = magnetic;
    if (!woundMagnetic.get_coil().get_turns_description()) {
        woundMagnetic.get_mutable_coil().wind();
    }

    // The same losses evaluated at ambient, independently of the loss loop.
    OperatingPoint operatingPointForWinding = inputs.get_operating_points()[0];
    double windingLossesAtAmbient = MagneticSimulator().calculate_winding_losses(operatingPointForWinding, woundMagnetic, ambientTemperature).get_winding_losses();
    MagnetizingInductance magnetizingInductanceModel(to_string(settings.get_reluctance_model()));
    magnetizingInductanceModel.set_topology(inputs.get_design_requirements().get_topology());
    OperatingPoint operatingPointForCore = inputs.get_operating_points()[0];
    auto fluxDensity = magnetizingInductanceModel.calculate_inductance_and_magnetic_flux_density(woundMagnetic, &operatingPointForCore).second;
    auto excitation = operatingPointForCore.get_excitations_per_winding()[0];
    excitation.set_magnetic_flux_density(fluxDensity);
    double coreLossesAtAmbient = CoreLosses().calculate_core_losses(woundMagnetic.get_core(), excitation, ambientTemperature).get_core_losses();

    MagneticSimulator simulator;
    simulator.set_topology(inputs.get_design_requirements().get_topology());
    OperatingPoint operatingPoint = inputs.get_operating_points()[0];

    SECTION("the standalone call, as the bindings make it: the core losses at ambient") {
        CoreLossesOutput coreLosses;
        REQUIRE_NOTHROW(coreLosses = simulator.calculate_core_losses(operatingPoint, woundMagnetic));
        INFO("core " << coreLosses.get_core_losses() << " W (at ambient " << coreLossesAtAmbient << " W) at " << coreLosses.get_temperature() << " C");
        CHECK(coreLosses.get_temperature() == ambientTemperature);
        CHECK_THAT(coreLosses.get_core_losses(), Catch::Matchers::WithinRel(coreLossesAtAmbient, 1e-12));
    }

    SECTION("asked for them, no temperature and the winding losses at ambient") {
        std::optional<WindingLossesOutput> windingLossesOutput;
        // Engaged beforehand: the call must empty it, not leave a stale temperature behind.
        std::optional<TemperatureOutput> temperatureOutput = TemperatureOutput();
        CoreLossesOutput coreLosses;
        REQUIRE_NOTHROW(coreLosses = simulator.calculate_core_losses(operatingPoint, woundMagnetic, &windingLossesOutput, &temperatureOutput));
        CHECK_FALSE(temperatureOutput);
        CHECK(coreLosses.get_temperature() == ambientTemperature);
        CHECK_THAT(coreLosses.get_core_losses(), Catch::Matchers::WithinRel(coreLossesAtAmbient, 1e-12));
        REQUIRE(windingLossesOutput);
        CHECK_THAT(windingLossesOutput->get_winding_losses(), Catch::Matchers::WithinRel(windingLossesAtAmbient, 1e-12));
    }

    SECTION("the advisers' policy still throws NoThermalOperatingPointException") {
        CHECK_THROWS_AS(simulator.calculate_core_losses(operatingPoint, woundMagnetic, nullptr, nullptr, WindingLossesScope::FULL,
                                                        NoThermalOperatingPointPolicy::THROW),
                        NoThermalOperatingPointException);
    }
    settings.reset();
}

TEST_CASE("simulate's JSON has no temperature key, not a null, for a point with no thermal operating point",
          "[processor][magnetic-simulator][abt-1876]") {
    // Alf decision 2026-10-14, ABT #1876: objects are schema-valid at every stage. On 28d80298 outputs[i].temperature
    // was written as an explicit null for the point with no thermal operating point.
    settings.reset();
    auto magnetic = load_web_default_design().first;
    json coolPoint = scaled_operating_point(0.01, "A hundredth of the current");
    json hotPoint = scaled_operating_point(1.0, "The web default");
    auto twoPoints = web_default_inputs_with(json::array({coolPoint, hotPoint}), magnetic);

    OpenMagnetics::Mas simulated;
    REQUIRE_NOTHROW(simulated = MagneticSimulator().simulate(twoPoints, magnetic));
    json masJson;
    to_json(masJson, simulated);
    REQUIRE(masJson["outputs"].size() == 2);
    REQUIRE(masJson["outputs"][0].contains("temperature"));
    CHECK(masJson["outputs"][0]["temperature"].is_object());
    CHECK_FALSE(masJson["outputs"][1].contains("temperature"));
    for (size_t index = 0; index < 2; ++index) {
        for (const auto& [key, value] : masJson["outputs"][index].items()) {
            INFO("outputs[" << index << "]." << key);
            CHECK_FALSE(value.is_null());
        }
    }

    // Read back as the web reads it: still no temperature at point 1, and analyze_results still names it.
    OpenMagnetics::Mas reread;
    from_json(masJson, reread);
    REQUIRE(reread.get_outputs().size() == 2);
    CHECK(reread.get_outputs()[0].get_temperature());
    CHECK_FALSE(reread.get_outputs()[1].get_temperature());
    auto findings = ResultsAnalyzer::analyze_results(reread);
    INFO("findings:" << describe_findings(findings));
    size_t noThermalAt1 = 0;
    for (const auto& finding : findings) {
        if (finding.code == ResultFindingCode::NO_THERMAL_OPERATING_POINT) {
            CHECK(finding.operatingPointIndex == 1);
            ++noThermalAt1;
        }
    }
    CHECK(noThermalAt1 == 1);
    settings.reset();
}

TEST_CASE("calculate_core_losses of a coil with turns and sections but no layers, as the web's builder calls it",
          "[processor][magnetic-simulator][abt-1886]") {
    // ABT #1886: the Magnetic Builder's call straight after Wire Advise (WebFrontend replay, models IGSE / Maniktala /
    // Zhang): one litz winding, 10 turns, sectionsDescription and turnsDescription but no layersDescription (the builder
    // keeps the advised coil's sections and turns, not its layers). It threw "[COIL_NOT_PROCESSED] ... Not wound by
    // layers" from WindingLosses::combine_turn_losses, and past that "Missing layers description" from the thermal
    // network's turn-to-turn insulation. The layers are now rebuilt from the sections before the solve.
    // Prod 426ec505 (core-only network): 0.00657 W, 25.8 C.
    settings.reset();
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "builder_core_losses_turns_without_layers.json");
    std::ifstream file(path);
    REQUIRE(file.is_open());
    json call = json::parse(file);
    json coreJson = call["magnetic"]["core"];
    // As WebLibMKF's calculate_core_losses builds its inputs.
    OpenMagnetics::Core core;
    core.set_functional_description(MAS::CoreFunctionalDescription(coreJson["functionalDescription"]));
    core.set_processed_description(MAS::CoreProcessedDescription(coreJson["processedDescription"]));
    core.set_geometrical_description(std::vector<MAS::CoreGeometricalDescriptionElement>(coreJson["geometricalDescription"]));
    OpenMagnetics::Coil coil(call["magnetic"]["coil"], false);
    REQUIRE(coil.get_turns_description());
    REQUIRE(coil.get_sections_description());
    REQUIRE_FALSE(coil.get_layers_description());
    double magnetizingInductance = MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(core, coil)
                                       .get_magnetizing_inductance().get_nominal().value();
    OpenMagnetics::Inputs inputs(call["inputs"], true, magnetizingInductance);
    auto operatingPoint = inputs.get_operating_point(call["operatingPointIndex"].get<size_t>());
    auto excitation = operatingPoint.get_excitations_per_winding()[0];
    if (!excitation.get_current()) {
        excitation.set_current(OpenMagnetics::Inputs::calculate_magnetizing_current(excitation, magnetizingInductance, true, 0.0));
        operatingPoint.get_mutable_excitations_per_winding()[0] = excitation;
    }
    auto coreLossesModelName = OpenMagnetics::defaults.coreLossesModelDefault;
    OpenMagnetics::from_json(call["models"]["coreLosses"], coreLossesModelName);
    auto reluctanceModelName = OpenMagnetics::defaults.reluctanceModelDefault;
    OpenMagnetics::from_json(call["models"]["gapReluctance"], reluctanceModelName);
    auto coreTemperatureModelName = OpenMagnetics::defaults.coreTemperatureModelDefault;
    OpenMagnetics::from_json(call["models"]["coreTemperature"], coreTemperatureModelName);
    settings.set_core_temperature_model(coreTemperatureModelName);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    MagneticSimulator simulator;
    simulator.set_core_losses_model_name(coreLossesModelName);
    simulator.set_reluctance_model_name(reluctanceModelName);

    std::optional<WindingLossesOutput> windingLossesOutput;
    std::optional<TemperatureOutput> temperatureOutput;
    CoreLossesOutput coreLosses;
    REQUIRE_NOTHROW(coreLosses = simulator.calculate_core_losses(operatingPoint, magnetic, &windingLossesOutput, &temperatureOutput));
    INFO("core " << coreLosses.get_core_losses() << " W at " << coreLosses.get_temperature() << " C (prod 426ec505: 0.00657 W, 25.8 C)");
    CHECK(std::isfinite(coreLosses.get_core_losses()));
    CHECK(coreLosses.get_core_losses() > 0);
    REQUIRE(windingLossesOutput);
    if (temperatureOutput) {
        INFO("hot spot " << temperatureOutput->get_maximum_temperature() << " C");
        CHECK(std::isfinite(temperatureOutput->get_maximum_temperature()));
        CHECK(temperatureOutput->get_maximum_temperature() > operatingPoint.get_conditions().get_ambient_temperature());
    }
    // Every breakdown, the per-layer one included: the layers the turns name were rebuilt.
    REQUIRE(windingLossesOutput->get_winding_losses_per_turn());
    CHECK(windingLossesOutput->get_winding_losses_per_turn()->size() == 10);
    REQUIRE(windingLossesOutput->get_winding_losses_per_section());
    CHECK(windingLossesOutput->get_winding_losses_per_section()->size() == 1);
    REQUIRE(windingLossesOutput->get_winding_losses_per_layer());
    REQUIRE(windingLossesOutput->get_winding_losses_per_layer()->size() == 2);
    CHECK(windingLossesOutput->get_winding_losses_per_layer()->at(0).get_name() == "Primary section 0 layer 0");
    CHECK(windingLossesOutput->get_winding_losses_per_layer()->at(1).get_name() == "Primary section 0 layer 1");
    // The caller's magnetic is left as given.
    CHECK_FALSE(magnetic.get_coil().get_layers_description());
    settings.reset();
}

// ABT #1879: three web-stored designs whose winding currents are 128 points on a CLOSED uniform
// axis ([0, T], last point = first). Their harmonics were taken as 128 OPEN samples (period one
// step long), so a pure sine carried fake content at every harmonic, and CoilMesher's phase check
// refused the design: "does not contain the harmonic its MAS harmonics list at 12600000 Hz" (web
// path) or 12800000 Hz (after mas_autocomplete). Both paths must simulate.
TEST_CASE("Web designs with 128-point closed-axis currents simulate", "[processor][magnetic-simulator][abt-1879]") {
    settings.reset();
    for (std::string name : {"simplest_transformer.json", "low_filling_transformer.json", "mas_autocomplete_bug.json"}) {
        std::ifstream file(OpenMagneticsTesting::get_test_data_path(std::source_location::current(), name));
        REQUIRE(file.good());
        auto masJson = nlohmann::json::parse(file);

        SECTION("web simulate: Inputs and Magnetic from json " + name) {
            OpenMagnetics::Inputs inputs(masJson["inputs"]);
            OpenMagnetics::Magnetic magnetic(masJson["magnetic"]);
            MagneticSimulator magneticSimulator;
            OpenMagnetics::Mas simulated;
            REQUIRE_NOTHROW(simulated = magneticSimulator.simulate(inputs, magnetic));
            CHECK(simulated.get_outputs().size() == inputs.get_operating_points().size());
        }
        SECTION("mas_autocomplete then simulate " + name) {
            OpenMagnetics::Mas mas(masJson);
            auto autocompleted = OpenMagnetics::mas_autocomplete(mas, false);
            nlohmann::json autocompletedJson;
            to_json(autocompletedJson, autocompleted);
            OpenMagnetics::Inputs inputs(autocompletedJson["inputs"]);
            OpenMagnetics::Magnetic magnetic(autocompletedJson["magnetic"]);
            MagneticSimulator magneticSimulator;
            OpenMagnetics::Mas simulated;
            REQUIRE_NOTHROW(simulated = magneticSimulator.simulate(inputs, magnetic));
            CHECK(simulated.get_outputs().size() == inputs.get_operating_points().size());
        }
    }
    settings.reset();
}

TEST_CASE("A coil with no primary winding is refused by name, not priced to NaN", "[processor][magnetic-simulator][abt-1843]") {
    // ABT #1843: one winding on isolationSide "secondary" and no primary. The winding losses refuse it
    // (no reference winding to refer the currents to); the core losses ahead of them used to end in
    // iGSE NaN on a negative flyback peak-to-peak instead.
    settings.reset();
    std::ifstream file(OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "low_filling_inductor.json"));
    REQUIRE(file.good());
    auto masJson = nlohmann::json::parse(file);
    OpenMagnetics::Inputs inputs(masJson["inputs"]);
    OpenMagnetics::Magnetic magnetic(masJson["magnetic"]);
    MagneticSimulator magneticSimulator;
    REQUIRE_THROWS_WITH(magneticSimulator.simulate(inputs, magnetic),
                        Catch::Matchers::ContainsSubstring("INVALID_COIL_CONFIGURATION") &&
                        Catch::Matchers::ContainsSubstring("no winding with isolationSide \"primary\""));
    settings.reset();
}

namespace {
    // ABT #1840 fixture: ETD 29/16/10 N87 with residual gaps, 10:10 turns of Round 0.5, 25 C, 100 kHz sinusoidal
    // voltage at `level` of B_sat(25 C), 1 A peak-to-peak.
    std::pair<OpenMagnetics::Magnetic, OpenMagnetics::Inputs> abt1840_design(double level) {
        const double frequency = 100e3;
        json magneticJson = {
            {"core", {{"name", "abt1840"}, {"functionalDescription", {{"type", "two-piece set"}, {"material", "N87"}, {"shape", "ETD 29/16/10"},
                      {"gapping", json::array({{{"type", "residual"}, {"length", 1e-5}}, {{"type", "residual"}, {"length", 1e-5}}, {{"type", "residual"}, {"length", 1e-5}}})},
                      {"numberStacks", 1}}}}},
            {"coil", {{"bobbin", "Basic"}, {"functionalDescription", json::array({
                      {{"name", "Primary"}, {"numberTurns", 10}, {"numberParallels", 1}, {"isolationSide", "primary"}, {"wire", "Round 0.5 - Grade 1"}},
                      {{"name", "Secondary"}, {"numberTurns", 10}, {"numberParallels", 1}, {"isolationSide", "secondary"}, {"wire", "Round 0.5 - Grade 1"}}})}}}};
        OpenMagnetics::Magnetic magnetic = magnetic_autocomplete(OpenMagnetics::Magnetic(magneticJson));
        auto material = magnetic.get_core().resolve_material();
        const double saturationFluxDensity = Core::get_magnetic_flux_density_saturation(material, 25.0, false);
        const double effectiveArea = magnetic.get_core().get_processed_description()->get_effective_parameters().get_effective_area();
        const double peakVoltage = 2 * M_PI * frequency * 10 * effectiveArea * level * saturationFluxDensity;
        auto excitation = [&](const std::string& name) {
            return json{{"name", name}, {"frequency", frequency},
                        {"voltage", {{"processed", {{"label", "sinusoidal"}, {"peakToPeak", 2 * peakVoltage}, {"offset", 0}, {"dutyCycle", 0.5}}}}},
                        {"current", {{"processed", {{"label", "sinusoidal"}, {"peakToPeak", 1.0}, {"offset", 0}, {"dutyCycle", 0.5}}}}}};
        };
        json inputsJson = {{"designRequirements", {{"magnetizingInductance", {{"nominal", 1e-3}}}, {"turnsRatios", json::array({{{"nominal", 1}}})}}},
                           {"operatingPoints", json::array({{{"name", "op"}, {"conditions", {{"ambientTemperature", 25.0}}},
                                                             {"excitationsPerWinding", json::array({excitation("Primary"), excitation("Secondary")})}}})}};
        return {magnetic, OpenMagnetics::Inputs(inputsJson)};
    }
}

TEST_CASE("Loss-temperature loop settles at the lowest stable operating point, not past the unstable one", "[processor][magnetic-simulator][abt-1840]") {
    // ABT #1840: the fixture above at 0.44 of B_sat(25 C). Its map g(T) (the core temperature the network returns
    // for the losses at T) has a stable fixed point near 97 C and an unstable one near 172 C; g(25 C) = 186 C. The
    // first step from ambient overshot the unstable point and the loop climbed to the data ceiling (208 C, 13 W
    // instead of 3 W) and held the core there with a data-limit finding. A part heated from ambient stops at the
    // lowest stable point.
    settings.reset();
    auto [magnetic, inputs] = abt1840_design(0.44);

    MagneticSimulator simulator;
    OpenMagnetics::Mas simulated;
    REQUIRE_NOTHROW(simulated = simulator.simulate(inputs, magnetic));
    REQUIRE(simulated.get_outputs().size() == 1);
    REQUIRE(simulated.get_outputs()[0].get_core_losses());
    const double coreTemperature = simulated.get_outputs()[0].get_core_losses()->get_temperature();
    INFO("core temperature " << coreTemperature << " C, core losses " << simulated.get_outputs()[0].get_core_losses()->get_core_losses() << " W");
    CHECK_THAT(coreTemperature, Catch::Matchers::WithinAbs(97, 5));
    // It got there through the search: the loop reached the ceiling first.
    CHECK(simulator.get_last_lowest_fixed_point_search_evaluations() > 0);
    for (const auto& finding : ResultsAnalyzer::analyze_results(simulated)) {
        CHECK(finding.code != ResultFindingCode::CORE_ABOVE_MATERIAL_DATA_LIMIT);
    }
    settings.reset();
}

TEST_CASE("A loss-temperature loop with no fixed point below the data limit pays only the search's coarse pass", "[processor][magnetic-simulator][abt-1840]") {
    // ABT #1840: the fixture above at 0.60 of B_sat(25 C), where h(T) = g(T) - T stays positive from
    // ambient to the data ceiling and the core is held there. The search for the lowest stable fixed point walks the
    // bracket from ambient to the loop's iterate next to the ceiling in ceil(sqrt(width / tolerance)) steps and evaluates the n points they reach.
    // The tolerance-spaced scan it replaced evaluated every 2 K of the bracket (63 networks here).
    settings.reset();
    auto [magnetic, inputs] = abt1840_design(0.60);

    MagneticSimulator simulator;
    OpenMagnetics::Mas simulated;
    REQUIRE_NOTHROW(simulated = simulator.simulate(inputs, magnetic));
    REQUIRE(simulated.get_outputs().size() == 1);
    REQUIRE(simulated.get_outputs()[0].get_core_losses());
    const double coreTemperature = simulated.get_outputs()[0].get_core_losses()->get_temperature();
    const double ambientTemperature = 25.0;
    INFO("core temperature " << coreTemperature << " C, search evaluations " << simulator.get_last_lowest_fixed_point_search_evaluations());
    // Held at the ceiling, with the finding that says so.
    bool heldAtDataLimit = false;
    for (const auto& finding : ResultsAnalyzer::analyze_results(simulated)) {
        heldAtDataLimit = heldAtDataLimit || finding.code == ResultFindingCode::CORE_ABOVE_MATERIAL_DATA_LIMIT;
    }
    REQUIRE(heldAtDataLimit);
    // The search walks from ambient to the loop's iterate next to the ceiling, where the core is held: at most that
    // wide (a core held at the end of a mu_i(T) table is reported at the ceiling, past the iterate).
    const double bracketWidth = coreTemperature - ambientTemperature;
    const size_t coarseSteps = static_cast<size_t>(std::ceil(std::sqrt(bracketWidth / MagneticSimulator::kOperatingTemperatureTolerance)));
    CHECK(simulator.get_last_lowest_fixed_point_search_evaluations() > 0);
    CHECK(simulator.get_last_lowest_fixed_point_search_evaluations() <= coarseSteps);
    settings.reset();
}

TEST_CASE("A loss-temperature loop that does not settle is that point's NO_THERMAL_OPERATING_POINT, not a refusal of the design",
          "[processor][magnetic-simulator][abt-1960]") {
    // Fair-Rite 78 EL 25/6.3, one winding of 6 turns on the primary side, 200 kHz, 25 C, completed as the web's
    // loader completes it (mas_autocomplete, not simulated). simulate threw for the whole design ("did not settle
    // within 100 iterations"); like a network that does not converge or a winding that runs away, it is that point
    // having no thermal operating point: losses at ambient, no temperature, and the finding from analyze_results.
    // The fixture was first that design as the web loaded it (single_winding_78_el25_on_primary.json): its core
    // losses jumped 4.19 -> 3.81 W where B(T) stepped down at 135.89 C. The unipolar DCM flux density read on the
    // magnetisation curve made that design settle (its core at 198.03 C, under the 200 C Curie point), so the path
    // was no longer reached. This fixture is the same design with its current scaled by 1.6 (3.6 A average, 9.18 A
    // peak; voltage unchanged). Its core heats towards the Curie point, where H_sat(T) goes to zero, and near
    // 197.97 C its losses step between 23.8 W and 0.11 W within a microkelvin, at the same peak flux density: the
    // core temperatures they produce, 387 C and 93 C, lie on either side of the step, so no temperature is within
    // the 2 K tolerance of the one its own losses produce and the damped iterates close in on the step from both
    // sides. If a change to the magnetizing or loss model removes that step, this fixture no longer reaches the
    // path and needs another.
    settings.reset();
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "single_winding_78_el25_on_primary_curie_step.json");
    auto completed = OpenMagnetics::mas_autocomplete(OpenMagneticsTesting::mas_loader(path), false, json{});
    auto magnetic = completed.get_magnetic();
    auto inputs = completed.get_inputs();
    REQUIRE(inputs.get_operating_points().size() == 1);
    const double ambientTemperature = inputs.get_operating_points()[0].get_conditions().get_ambient_temperature();
    REQUIRE(ambientTemperature == 25);

    SECTION("the advisers' policy still throws, as no thermal operating point") {
        OpenMagnetics::Magnetic woundMagnetic = magnetic;
        REQUIRE(woundMagnetic.get_coil().get_turns_description());
        OperatingPoint operatingPoint = inputs.get_operating_points()[0];
        std::optional<WindingLossesOutput> windingLossesOutput;
        std::optional<TemperatureOutput> temperatureOutput;
        CHECK_THROWS_MATCHES(MagneticSimulator().calculate_core_losses(operatingPoint, woundMagnetic, &windingLossesOutput, &temperatureOutput,
                                                                       WindingLossesScope::FULL, NoThermalOperatingPointPolicy::THROW),
                             NoThermalOperatingPointException,
                             Catch::Matchers::MessageMatches(Catch::Matchers::ContainsSubstring("did not settle within")));
    }

    SECTION("simulate gives the point its losses at ambient, and analyze_results names it") {
        MagneticSimulator simulator;
        OpenMagnetics::Mas simulated;
        REQUIRE_NOTHROW(simulated = simulator.simulate(inputs, magnetic));
        REQUIRE(simulated.get_outputs().size() == 1);
        const auto& output = simulated.get_outputs()[0];
        CHECK_FALSE(output.get_temperature());
        REQUIRE(output.get_core_losses());
        REQUIRE(output.get_winding_losses());
        CHECK(output.get_core_losses()->get_temperature() == ambientTemperature);

        auto findings = ResultsAnalyzer::analyze_results(simulated);
        INFO("findings:" << describe_findings(findings));
        std::vector<ResultFinding> noThermal;
        for (const auto& finding : findings) {
            if (finding.code == ResultFindingCode::NO_THERMAL_OPERATING_POINT) {
                noThermal.push_back(finding);
            }
        }
        REQUIRE(noThermal.size() == 1);
        CHECK(noThermal[0].operatingPointIndex == 0);
        CHECK_THAT(noThermal[0].message, Catch::Matchers::ContainsSubstring("No steady temperature") &&
                                         Catch::Matchers::ContainsSubstring("did not settle within") &&
                                         Catch::Matchers::ContainsSubstring("at the ambient temperature"));
        REQUIRE(noThermal[0].details);
        CHECK_THAT(noThermal[0].details.value(), Catch::Matchers::ContainsSubstring("History:"));
    }
    settings.reset();
}


namespace {
    // Settings::simulatorLossTemperatureLoop is global: every test that turns it off restores it on every exit
    // path, a failed REQUIRE included.
    struct LossTemperatureLoopOff {
        SettingsGuard<bool> guard{Settings::GetInstance(), &Settings::get_simulator_loss_temperature_loop,
                                  &Settings::set_simulator_loss_temperature_loop, false};
        ~LossTemperatureLoopOff() { Settings::GetInstance().reset(); }
    };

    std::string simulated_json(const OpenMagnetics::Inputs& inputs, const OpenMagnetics::Magnetic& magnetic) {
        json masJson;
        to_json(masJson, MagneticSimulator().simulate(inputs, magnetic));
        return masJson.dump();
    }
}

// Alf, 2026-10-22: the loss-temperature loop is a setting (Settings::simulatorLossTemperatureLoop), on by default.
// On, simulate is what it was before the setting existed: setting it on explicitly changes no byte of its output.
TEST_CASE("The loss-temperature loop is on by default, and on is the simulate it always was",
          "[processor][magnetic-simulator][abt-2003]") {
    settings.reset();
    REQUIRE(settings.get_simulator_loss_temperature_loop());
    auto path = get_examples_dir() / "01_simple_inductor_etd34_n87.json";
    REQUIRE(std::filesystem::exists(path));
    auto mas = OpenMagneticsTesting::mas_loader(path.string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
    const std::string byDefault = simulated_json(inputs, magnetic);
    std::string explicitlyOn;
    {
        SettingsGuard<bool> on(settings, &Settings::get_simulator_loss_temperature_loop, &Settings::set_simulator_loss_temperature_loop, true);
        explicitlyOn = simulated_json(inputs, magnetic);
    }
    CHECK(byDefault == explicitlyOn);
    {
        LossTemperatureLoopOff off;
        CHECK_FALSE(settings.get_simulator_loss_temperature_loop());
        CHECK(simulated_json(inputs, magnetic) != byDefault);
    }
    // reset() turns it back on.
    CHECK(settings.get_simulator_loss_temperature_loop());
    settings.set_simulator_loss_temperature_loop(false);
    settings.reset();
    CHECK(settings.get_simulator_loss_temperature_loop());
}

// Alf, 2026-10-22: with the loop off, simulate evaluates every loss at the operating point's temperature and still
// publishes a temperature: the hottest point the thermal network gives ONCE for exactly those losses, fed back
// into nothing. There is no loop, so no finding about one.
TEST_CASE("With the loss-temperature loop off, simulate gives the losses at the operating temperature and one thermal pass on them",
          "[processor][magnetic-simulator][abt-2003]") {
    settings.reset();
    const double tolerance = MagneticSimulator::kOperatingTemperatureTolerance;
    auto path = get_examples_dir() / "01_simple_inductor_etd34_n87.json";
    REQUIRE(std::filesystem::exists(path));
    auto mas = OpenMagneticsTesting::mas_loader(path.string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
    if (!magnetic.get_coil().get_turns_description()) {
        magnetic.get_mutable_coil().wind();
    }
    const auto looped = MagneticSimulator().simulate(inputs, magnetic);

    LossTemperatureLoopOff off;
    MagneticSimulator simulator;
    const auto simulated = simulator.simulate(inputs, magnetic);
    REQUIRE(simulated.get_outputs().size() == inputs.get_operating_points().size());
    for (size_t index = 0; index < simulated.get_outputs().size(); ++index) {
        INFO("operating point " << index);
        const auto& output = simulated.get_outputs()[index];
        const auto& operatingPoint = simulated.get_inputs().get_operating_points()[index];
        const double ambientTemperature = inputs.get_operating_points()[index].get_conditions().get_ambient_temperature();
        REQUIRE(operatingPoint.get_conditions().get_ambient_temperature() == ambientTemperature);

        // The losses: at the operating temperature, as the models give them there.
        REQUIRE(output.get_core_losses());
        REQUIRE(output.get_winding_losses());
        CHECK(output.get_core_losses()->get_temperature() == ambientTemperature);
        REQUIRE(output.get_winding_losses()->get_temperature());
        CHECK(output.get_winding_losses()->get_temperature().value() == ambientTemperature);
        auto excitation = operatingPoint.get_excitations_per_winding()[0];
        REQUIRE(excitation.get_magnetic_flux_density());
        const double coreLossesAtAmbient = CoreLosses().calculate_core_losses(magnetic.get_core(), excitation, ambientTemperature).get_core_losses();
        OperatingPoint operatingPointForWinding = inputs.get_operating_points()[index];
        const double windingLossesAtAmbient =
            MagneticSimulator().calculate_winding_losses(operatingPointForWinding, magnetic, ambientTemperature).get_winding_losses();
        CAPTURE(output.get_core_losses()->get_core_losses(), coreLossesAtAmbient, output.get_winding_losses()->get_winding_losses(),
                windingLossesAtAmbient);
        CHECK_THAT(output.get_core_losses()->get_core_losses(), Catch::Matchers::WithinRel(coreLossesAtAmbient, 1e-12));
        CHECK_THAT(output.get_winding_losses()->get_winding_losses(), Catch::Matchers::WithinRel(windingLossesAtAmbient, 1e-12));

        // The temperature: one pass of the thermal network on exactly those losses.
        REQUIRE(output.get_temperature());
        auto onePass = Temperature(magnetic, TemperatureConfig::fromSimulatedOutput(operatingPoint, output)).calculateTemperatures();
        CAPTURE(output.get_temperature()->get_maximum_temperature(), onePass.maximumTemperature);
        CHECK_THAT(output.get_temperature()->get_maximum_temperature(), Catch::Matchers::WithinRel(onePass.maximumTemperature, 1e-12));
        CHECK(output.get_temperature()->get_initial_temperature() == ambientTemperature);
        // A real rise, so the pass is measured: the losses heat the part well past the loop's tolerance ...
        CHECK(onePass.maximumTemperature > ambientTemperature + tolerance);
        // ... and they are not the loop's: the loop evaluates them at that higher temperature.
        const auto& loopedOutput = looped.get_outputs()[index];
        CAPTURE(loopedOutput.get_core_losses()->get_temperature(), loopedOutput.get_core_losses()->get_core_losses());
        CHECK(loopedOutput.get_core_losses()->get_temperature() > ambientTemperature + tolerance);
        CHECK(std::fabs(loopedOutput.get_core_losses()->get_core_losses() / output.get_core_losses()->get_core_losses() - 1) > 1e-3);
    }

    // No loop, no finding about one.
    auto findings = ResultsAnalyzer::analyze_results(simulated);
    INFO("findings:" << describe_findings(findings));
    for (const auto& finding : findings) {
        CHECK(finding.code != ResultFindingCode::NO_THERMAL_OPERATING_POINT);
        CHECK(finding.code != ResultFindingCode::CORE_ABOVE_MATERIAL_DATA_LIMIT);
    }
}

// The web's default design has no thermal operating point with the loop on (NO_THERMAL_OPERATING_POINT, ABT #1876):
// the loop runs away from its first iterate. With the loop off that first iterate is the answer: losses at ambient
// and the network's hottest point for them, published, and no NO_THERMAL_OPERATING_POINT finding.
TEST_CASE("With the loss-temperature loop off, the web's default design gets its one-pass temperature and no loop finding",
          "[processor][magnetic-simulator][abt-2003]") {
    settings.reset();
    auto [magnetic, inputs] = load_web_default_design();
    REQUIRE(inputs.get_operating_points().size() == 1);
    const double ambientTemperature = inputs.get_operating_points()[0].get_conditions().get_ambient_temperature();
    const auto looped = MagneticSimulator().simulate(inputs, magnetic);
    REQUIRE(looped.get_outputs().size() == 1);
    // The loop's answer for this design: no temperature (see the ABT #1876 test above).
    REQUIRE_FALSE(looped.get_outputs()[0].get_temperature());

    LossTemperatureLoopOff off;
    OpenMagnetics::Mas simulated;
    REQUIRE_NOTHROW(simulated = MagneticSimulator().simulate(inputs, magnetic));
    REQUIRE(simulated.get_outputs().size() == 1);
    const auto& output = simulated.get_outputs()[0];
    // Same losses as the loop's no-thermal-operating-point answer: both are the losses at ambient.
    CHECK(output.get_core_losses()->get_temperature() == ambientTemperature);
    CHECK_THAT(output.get_core_losses()->get_core_losses(),
               Catch::Matchers::WithinRel(looped.get_outputs()[0].get_core_losses()->get_core_losses(), 1e-12));
    CHECK_THAT(output.get_winding_losses()->get_winding_losses(),
               Catch::Matchers::WithinRel(looped.get_outputs()[0].get_winding_losses()->get_winding_losses(), 1e-12));
    REQUIRE(output.get_temperature());
    OpenMagnetics::Magnetic woundMagnetic = magnetic;
    if (!woundMagnetic.get_coil().get_turns_description()) {
        woundMagnetic.get_mutable_coil().wind();
    }
    auto onePass = Temperature(woundMagnetic, TemperatureConfig::fromSimulatedOutput(simulated.get_inputs().get_operating_points()[0], output))
                       .calculateTemperatures();
    CAPTURE(output.get_temperature()->get_maximum_temperature(), onePass.maximumTemperature);
    CHECK_THAT(output.get_temperature()->get_maximum_temperature(), Catch::Matchers::WithinRel(onePass.maximumTemperature, 1e-12));

    auto findings = ResultsAnalyzer::analyze_results(simulated);
    INFO("findings:" << describe_findings(findings));
    for (const auto& finding : findings) {
        CHECK(finding.code != ResultFindingCode::NO_THERMAL_OPERATING_POINT);
        CHECK(finding.code != ResultFindingCode::CORE_ABOVE_MATERIAL_DATA_LIMIT);
    }
}

// Advisers keep their own behaviour: the setting is for simulate and the builder. While an adviser exists the
// effective setting reads on (Settings::get_simulator_loss_temperature_loop), so the loss solves its filters and
// its simulate run are the loop's, bit for bit.
TEST_CASE("Under an adviser the loss solve runs the loss-temperature loop whatever the setting says",
          "[processor][magnetic-simulator][abt-2003]") {
    settings.reset();
    auto path = get_examples_dir() / "01_simple_inductor_etd34_n87.json";
    REQUIRE(std::filesystem::exists(path));
    auto mas = OpenMagneticsTesting::mas_loader(path.string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
    if (!magnetic.get_coil().get_turns_description()) {
        magnetic.get_mutable_coil().wind();
    }
    const double ambientTemperature = inputs.get_operating_points()[0].get_conditions().get_ambient_temperature();
    OperatingPoint loopedOperatingPoint = inputs.get_operating_points()[0];
    const auto looped = MagneticSimulator().calculate_core_losses(loopedOperatingPoint, magnetic);
    REQUIRE(looped.get_temperature() > ambientTemperature + MagneticSimulator::kOperatingTemperatureTolerance);

    LossTemperatureLoopOff off;
    {
        MagneticAdviser adviser;
        CHECK(settings.get_simulator_loss_temperature_loop());
        OperatingPoint operatingPoint = inputs.get_operating_points()[0];
        const auto underAdviser = MagneticSimulator().calculate_core_losses(operatingPoint, magnetic);
        CHECK(underAdviser.get_temperature() == looped.get_temperature());
        CHECK(underAdviser.get_core_losses() == looped.get_core_losses());
    }
    // The adviser gone, the caller's choice reads back.
    CHECK_FALSE(settings.get_simulator_loss_temperature_loop());
    OperatingPoint operatingPoint = inputs.get_operating_points()[0];
    CHECK(MagneticSimulator().calculate_core_losses(operatingPoint, magnetic).get_temperature() == ambientTemperature);
}


namespace {
void strip_nulls_abt1942_simulator(json& value) {
    if (value.is_array()) {
        for (auto& element : value) {
            strip_nulls_abt1942_simulator(element);
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
            strip_nulls_abt1942_simulator(*it);
            ++it;
        }
    }
}
}

TEST_CASE("Loss loop whose heated core cannot carry the voltage's flux on its curve has no thermal operating point",
          "[processor][magnetic-simulator][abt-1942]") {
    // ABT #1942 (Alf, 2026-10-19): a voltage-driven DCM point has its flux from Faraday and its magnetizing current read
    // back from the material's magnetisation curve at the core's temperature. slow_simulation's second operating point
    // (flyback, material 95, Curie 220 C) is fine at 25 C ambient (core ~110 C), but from a hotter ambient its losses
    // drive the core to where 95's curve cannot carry that flux: past its end above B_sat(T) (the readback's
    // SaturationException) or, closer to Curie where H_sat(T) goes to zero and the curve ends at 1000 H_sat(T), past
    // its end while still under B_sat(T) (MagnetisationCurveExhaustedException). Either refusal escaped the loss solve
    // and took the MagneticAdviser's whole call with it (error_paint_field, catalogue_missing_turns,
    // slow_simulation). At a heated iterate it is that point's NoThermalOperatingPointException (ABT #1876, #1960):
    // losses at ambient, no temperature, and analyze_results' NO_THERMAL_OPERATING_POINT.
    settings.reset();
    std::ifstream file(OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "slow_simulation.json"));
    REQUIRE(file.good());
    json masJson = json::parse(file);
    strip_nulls_abt1942_simulator(masJson);
    OpenMagnetics::Magnetic magnetic(masJson["magnetic"]);
    if (!magnetic.get_coil().get_turns_description()) {
        magnetic.get_mutable_coil().wind();
    }
    OpenMagnetics::Inputs allInputs(masJson["inputs"]);
    REQUIRE(allInputs.get_operating_points().size() > 1);

    // 85 C: past the curve's end above B_sat; 115 C: past its end under B_sat (near Curie).
    for (double ambientTemperature : {85.0, 115.0}) {
        DYNAMIC_SECTION("ambient " << ambientTemperature << " C") {
            OperatingPoint operatingPoint = allInputs.get_operating_points()[1];
            operatingPoint.get_mutable_conditions().set_ambient_temperature(ambientTemperature);
            const std::string refusal = ambientTemperature < 100 ? "beyond the end of the magnetisation curve" : "reaches only";
            {
                OperatingPoint thrown = operatingPoint;
                std::optional<WindingLossesOutput> windingLossesOutput;
                std::optional<TemperatureOutput> temperatureOutput;
                CHECK_THROWS_MATCHES(MagneticSimulator().calculate_core_losses(thrown, magnetic, &windingLossesOutput, &temperatureOutput,
                                                                               WindingLossesScope::OHMIC, NoThermalOperatingPointPolicy::THROW),
                                     NoThermalOperatingPointException,
                                     Catch::Matchers::MessageMatches(Catch::Matchers::ContainsSubstring("No steady temperature: heated by its own losses to") &&
                                                                     Catch::Matchers::ContainsSubstring(refusal)));
            }
            OpenMagnetics::Inputs inputs = allInputs;
            inputs.set_operating_points({operatingPoint});
            MagneticSimulator simulator;
            OpenMagnetics::Mas simulated;
            REQUIRE_NOTHROW(simulated = simulator.simulate(inputs, magnetic));
            REQUIRE(simulated.get_outputs().size() == 1);
            const auto& output = simulated.get_outputs()[0];
            CHECK_FALSE(output.get_temperature());
            REQUIRE(output.get_core_losses());
            CHECK(output.get_core_losses()->get_temperature() == ambientTemperature);
            auto findings = ResultsAnalyzer::analyze_results(simulated);
            INFO("findings:" << describe_findings(findings));
            size_t noThermal = 0;
            for (const auto& finding : findings) {
                if (finding.code == ResultFindingCode::NO_THERMAL_OPERATING_POINT) {
                    ++noThermal;
                    CHECK_THAT(finding.message, Catch::Matchers::ContainsSubstring(refusal));
                }
            }
            CHECK(noThermal == 1);
        }
    }
    settings.reset();
}


// ABT #2001: the winding-temperature finding is one plain sentence with rounded numbers. It read "the design's
// hottest point is 2474.39 C, above winding 0 ('Primary')'s the thermal class of its coating material
// Polyurethane 155 of 155 C".
TEST_CASE("Analyze_results states a winding above its coating's thermal class in one rounded sentence",
          "[processor][magnetic-simulator][abt-2001]") {
    settings.reset();
    auto path = std::filesystem::path{std::source_location::current().file_name()}
                    .parent_path().append("..").append("MAS").append("examples").append("02_flyback_efd25_3c95.json");
    REQUIRE(std::filesystem::exists(path));
    auto mas = OpenMagneticsTesting::mas_loader(path.string());
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    // The primary's enamel names its material, whose thermal class is 155 C.
    auto& primary = magnetic.get_mutable_coil().get_mutable_functional_description()[0];
    json wireJson;
    to_json(wireJson, OpenMagnetics::Coil::resolve_wire(primary));
    wireJson["coating"]["material"] = "Polyurethane 155";
    wireJson["coating"].erase("temperatureRating");
    primary.set_wire(OpenMagnetics::Wire(wireJson));

    MagneticSimulator magneticSimulator;
    OpenMagnetics::Mas simulated;
    REQUIRE_NOTHROW(simulated = magneticSimulator.simulate(mas.get_inputs(), magnetic));
    REQUIRE(simulated.get_outputs().size() == 1);
    auto output = simulated.get_outputs()[0];
    REQUIRE(output.get_temperature());
    auto temperature = output.get_temperature().value();
    temperature.set_maximum_temperature(2474.39);
    output.set_temperature(temperature);

    auto findings = ResultsAnalyzer::analyze_temperature(simulated.get_magnetic(), simulated.get_inputs().get_operating_points()[0], output, 0);
    std::vector<ResultFinding> windingFindings;
    std::string all;
    for (const auto& finding : findings) {
        all += "\n" + to_string(finding.code) + ": " + finding.message;
        if (finding.code == ResultFindingCode::WINDING_TEMPERATURE_ABOVE_LIMIT) {
            windingFindings.push_back(finding);
        }
    }
    INFO("findings:" << all);
    REQUIRE(!windingFindings.empty());
    CHECK_THAT(windingFindings[0].message, Catch::Matchers::EndsWith(
        ": the hottest point, 2474 C, is above the 155 C thermal class of the coating material (Polyurethane 155) of winding 'Primary'. "
        "Lower the losses, improve the cooling, or use a wire rated for a higher temperature."));
}

namespace {
    // ABT #2052: slow_simulation (flyback, material 95, Curie 220 C) whose second operating point ('Max. input volt.')
    // is simulated twice: at its own 25 C, where it is fine, and at an ambient so hot that its voltage's flux density,
    // read on 95's magnetisation curve at that ambient, lies past the curve's end: at 160 C above B_sat
    // (SaturationException: 0.404 T asked, the curve ends at 0.393 T), at 190 C under it, where the curve ends below
    // B_sat (MagnetisationCurveExhaustedException: 0.225 T asked, the curve ends at 0.219 T).
    struct PastCurveAtAmbientDesign {
        OpenMagnetics::Magnetic magnetic;
        OpenMagnetics::Inputs inputs;
    };
    PastCurveAtAmbientDesign past_curve_at_ambient_design(double hotAmbientTemperature) {
        std::ifstream file(OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "slow_simulation.json"));
        REQUIRE(file.good());
        json masJson = json::parse(file);
        strip_nulls_abt1942_simulator(masJson);
        PastCurveAtAmbientDesign design{OpenMagnetics::Magnetic(masJson["magnetic"]), OpenMagnetics::Inputs(masJson["inputs"])};
        if (!design.magnetic.get_coil().get_turns_description()) {
            design.magnetic.get_mutable_coil().wind();
        }
        // Its second operating point twice: as given (25 C, where its loss loop settles) and at the hot ambient.
        REQUIRE(design.inputs.get_operating_points().size() == 2);
        OperatingPoint hot = design.inputs.get_operating_points()[1];
        REQUIRE(hot.get_conditions().get_ambient_temperature() == 25);
        hot.get_mutable_conditions().set_ambient_temperature(hotAmbientTemperature);
        design.inputs.set_operating_points({design.inputs.get_operating_points()[1], hot});
        return design;
    }
    // The two ambients of the guards below and the refusal each meets at it.
    const std::vector<std::pair<double, std::string>> pastCurveAtAmbientCases = {
        {160.0, "beyond the end of the magnetisation curve"},
        {190.0, "reaches only"}};
}

TEST_CASE("A point whose flux density lies past the end of the magnetisation curve at ambient is that point's CORE_SATURATES",
          "[processor][magnetic-simulator][abt-2052]") {
    // ABT #2052 (Alf, 2026-10-24): since ABT #1942 the magnetizing model refuses a flux density past the end of the
    // material's curve; at a heated iterate of the loss solve that is the point's NoThermalOperatingPointException, but at
    // ambient it propagated and failed the whole simulation. It is that point's finding: simulate leaves its
    // flux-dependent outputs (inductance, core losses, temperature) unset, gives it its winding losses at ambient,
    // simulates every other point, and analyze_results reports CORE_SATURATES for it.
    settings.reset();
    for (const auto& [hotAmbientTemperature, refusal] : pastCurveAtAmbientCases) {
        DYNAMIC_SECTION("ambient " << hotAmbientTemperature << " C") {
            auto design = past_curve_at_ambient_design(hotAmbientTemperature);
            {
                // The precondition: at that ambient the magnetizing model cannot read the point.
                OperatingPoint operatingPoint = design.inputs.get_operating_points()[1];
                MagnetizingInductance magnetizingInductanceModel;
                magnetizingInductanceModel.set_topology(design.inputs.get_design_requirements().get_topology());
                CHECK_THROWS_WITH(magnetizingInductanceModel.calculate_inductance_and_magnetic_flux_density(design.magnetic, &operatingPoint),
                                  Catch::Matchers::ContainsSubstring(refusal));
            }
            MagneticSimulator simulator;
            OpenMagnetics::Mas simulated;
            REQUIRE_NOTHROW(simulated = simulator.simulate(design.inputs, design.magnetic));
            REQUIRE(simulated.get_outputs().size() == 2);

            // The other point is complete.
            const auto& complete = simulated.get_outputs()[0];
            REQUIRE(complete.get_inductance());
            REQUIRE(complete.get_core_losses());
            CHECK(complete.get_core_losses()->get_magnetic_flux_density());
            CHECK(complete.get_core_losses()->get_core_losses() > 0);
            CHECK(complete.get_temperature());
            REQUIRE(complete.get_winding_losses());
            CHECK(complete.get_winding_losses()->get_winding_losses() > 0);

            // The saturated point: its flux-dependent outputs unset, its winding losses at ambient.
            const auto& saturated = simulated.get_outputs()[1];
            CHECK_FALSE(saturated.get_inductance());
            CHECK_FALSE(saturated.get_core_losses());
            CHECK_FALSE(saturated.get_temperature());
            REQUIRE(saturated.get_winding_losses());
            CHECK(saturated.get_winding_losses()->get_winding_losses() > 0);
            REQUIRE(saturated.get_winding_losses()->get_temperature());
            CHECK(saturated.get_winding_losses()->get_temperature().value() == hotAmbientTemperature);
            // No magnetizing current was derived for it: its excitation is as given.
            const auto& givenExcitation = design.inputs.get_operating_points()[1].get_excitations_per_winding()[0];
            const auto& simulatedExcitation = simulated.get_inputs().get_operating_points()[1].get_excitations_per_winding()[0];
            CHECK(simulatedExcitation.get_magnetizing_current().has_value() == givenExcitation.get_magnetizing_current().has_value());
            CHECK_FALSE(simulatedExcitation.get_magnetic_flux_density());
            // Its temperature cannot be computed from that output.
            CHECK_THROWS_AS(TemperatureConfig::fromSimulatedOutput(simulated.get_inputs().get_operating_points()[1], saturated),
                            InvalidInputException);

            auto findings = ResultsAnalyzer::analyze_results(simulated);
            INFO("findings:" << describe_findings(findings));
            size_t coreSaturates = 0;
            for (const auto& finding : findings) {
                CHECK_THAT(finding.message, !Catch::Matchers::ContainsSubstring("ABT #"));
                if (finding.operatingPointIndex && finding.operatingPointIndex.value() == 1) {
                    CHECK(finding.code != ResultFindingCode::NO_THERMAL_OPERATING_POINT);
                }
                // Point 0 has its own CORE_SATURATES (its heated core, ~141 C, is past its knee); only point 1's is counted.
                if (finding.code == ResultFindingCode::CORE_SATURATES && finding.operatingPointIndex && finding.operatingPointIndex.value() == 1) {
                    ++coreSaturates;
                    CHECK_THAT(finding.message, Catch::Matchers::ContainsSubstring(refusal) &&
                                                Catch::Matchers::ContainsSubstring("are not calculated for this operating point"));
                }
            }
            CHECK(coreSaturates == 1);

            // An output without core losses at a point the model DOES read is not this case: analyze_results throws.
            auto edited = simulated;
            auto outputs = edited.get_outputs();
            outputs[0].set_core_losses(std::nullopt);
            edited.set_outputs(outputs);
            CHECK_THROWS_AS(ResultsAnalyzer::analyze_results(edited), InvalidInputException);
        }
    }
    settings.reset();
}

TEST_CASE("mas_autocomplete completes a design with a point past the end of the magnetisation curve at ambient",
          "[processor][magnetic-simulator][autocomplete][abt-2052]") {
    // ABT #2052: autocomplete never refuses a design for a saturating point. It fills the design, derives no magnetizing
    // current for that point, and its simulation leaves that point's flux-dependent outputs unset.
    settings.reset();
    for (const auto& [hotAmbientTemperature, refusal] : pastCurveAtAmbientCases) {
        DYNAMIC_SECTION("ambient " << hotAmbientTemperature << " C") {
            auto design = past_curve_at_ambient_design(hotAmbientTemperature);
            // As the web hands it over: the operating points without the magnetizing currents autocomplete derives.
            auto operatingPoints = design.inputs.get_operating_points();
            for (auto& operatingPoint : operatingPoints) {
                for (auto& excitation : operatingPoint.get_mutable_excitations_per_winding()) {
                    excitation.set_magnetizing_current(std::nullopt);
                }
            }
            design.inputs.set_operating_points(operatingPoints);
            OpenMagnetics::Mas mas;
            mas.set_inputs(design.inputs);
            mas.set_magnetic(design.magnetic);
            OpenMagnetics::Mas completed;
            REQUIRE_NOTHROW(completed = OpenMagnetics::mas_autocomplete(mas));
            REQUIRE(completed.get_outputs().size() == 2);
            CHECK(completed.get_outputs()[0].get_core_losses());
            CHECK(completed.get_outputs()[0].get_temperature());
            CHECK_FALSE(completed.get_outputs()[1].get_core_losses());
            CHECK_FALSE(completed.get_outputs()[1].get_temperature());
            CHECK(completed.get_outputs()[1].get_winding_losses());
            const auto& completedOperatingPoints = completed.get_inputs().get_operating_points();
            REQUIRE(completedOperatingPoints.size() == 2);
            CHECK(completedOperatingPoints[0].get_excitations_per_winding()[0].get_magnetizing_current());
            CHECK_FALSE(completedOperatingPoints[1].get_excitations_per_winding()[0].get_magnetizing_current());

            // Without simulating as well.
            REQUIRE_NOTHROW(OpenMagnetics::mas_autocomplete(mas, false));
        }
    }
    settings.reset();
}

// Alf, 2026-10-27: a powder has no knee at the field its saturation entry is tabulated at (MAS lists 294 powder
// entries at 100 Oe), so the peak-field criterion of check_saturation does not apply to it; a ferrite keeps it.
// CSC High Flux 60 is the material of the published Vienna 30 kW boost choke, which ran past 100 Oe.
TEST_CASE("A powder driven past the field its saturation is tabulated at is not saturated by that field alone", "[powder-knee]") {
    const double temperature = 25;
    auto powder = find_core_material_by_name("CSC High Flux 60");
    REQUIRE(powder.get_material() == MAS::MaterialType::POWDER);
    const double powderSaturationField = Core::get_magnetic_field_strength_saturation(powder, temperature);
    const double powderField = 2.5 * powderSaturationField;
    const double powderPeak = InitialPermeability::get_flux_density_for_magnetic_field_dc_bias(powder, powderField, temperature, 100e3);
    INFO("CSC High Flux 60: H " << powderField << " A/m (H_sat " << powderSaturationField << "), B " << powderPeak << " T");
    REQUIRE(powderPeak < Core::get_magnetic_flux_density_saturation(powder, temperature, false));
    DcBiasPoint powderPoint{powderField, temperature, 100e3, powderField};
    auto powderVerdict = ResultsAnalyzer::check_saturation(powder, powderPeak, temperature, 1.0, powderPoint);
    INFO(powderVerdict.message);
    REQUIRE(powderVerdict.biasedPermeability);
    CHECK(powderVerdict.biasedPermeability.value() > 1);
    CHECK_FALSE(powderVerdict.fieldStrengthAboveSaturation);
    CHECK_FALSE(powderVerdict.saturates);

    auto ferrite = find_core_material_by_name("3C95");
    const double ferriteSaturationField = Core::get_magnetic_field_strength_saturation(ferrite, temperature);
    DcBiasPoint ferritePoint{0, temperature, 100e3, 2.5 * ferriteSaturationField};
    auto ferriteVerdict = ResultsAnalyzer::check_saturation(ferrite, 0.3, temperature, 1.0, ferritePoint);
    INFO(ferriteVerdict.message);
    CHECK(ferriteVerdict.fieldStrengthAboveSaturation);
    CHECK(ferriteVerdict.saturates);
}

namespace {
// A two-piece set whose data states its gapping as catalogue exports write an ungapped part: a subtractive gap of
// length 0 on the centre column and the residual gaps of the outer columns (five Wurth catalogue parts carry it).
OpenMagnetics::Magnetic zero_length_subtractive_gap_magnetic(const std::string& coreType) {
    json gapping = json::array();
    gapping.push_back({{"type", "subtractive"}, {"length", 0.0}});
    gapping.push_back({{"type", "residual"}, {"length", 5e-6}});
    gapping.push_back({{"type", "residual"}, {"length", 5e-6}});
    json coreJson;
    coreJson["functionalDescription"]["name"] = "Zero length gap";
    coreJson["functionalDescription"]["type"] = coreType;
    coreJson["functionalDescription"]["material"] = "3C95";
    coreJson["functionalDescription"]["shape"] = "EFD 20/10/7";
    coreJson["functionalDescription"]["gapping"] = gapping;
    coreJson["functionalDescription"]["numberStacks"] = 1;
    OpenMagnetics::Core core(coreJson);
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, false);
    json bobbinJson;
    to_json(bobbinJson, bobbin);
    json coilJson;
    coilJson["bobbin"] = bobbinJson;
    coilJson["functionalDescription"] = json::array();
    json windingJson;
    windingJson["name"] = "Primary";
    windingJson["numberTurns"] = 12;
    windingJson["numberParallels"] = 1;
    windingJson["isolationSide"] = "primary";
    windingJson["wire"] = "Round 0.475 - Grade 1";
    coilJson["functionalDescription"].push_back(windingJson);
    OpenMagnetics::Coil coil(coilJson, 1, WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                             CoilAlignment::CENTERED, CoilAlignment::CENTERED);
    coil.delimit_and_compact();
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    return magnetic;
}
}  // namespace

TEST_CASE("Test_Simulate_Two_Piece_Set_With_Zero_Length_Subtractive_Gap_Is_Ungapped_With_A_Warning",
          "[simulator][gapping][zero-length-gap]") {
    // ABT #2124 (Alf, 2026-10-28): "a 0 gap length in a two piece set should be equivalent to ungapped/residual,
    // with just a warning, not an error". The neutral-window field (the default winding-loss field of a single
    // rectangular window) refused that gap in the final simulation ("a gap at y = 0 opens nowhere onto the winding
    // window"), so those parts dropped out of every catalogue search at its last step.
    settings.reset();
    Logger::getInstance().enableCollector(LogLevel::WARNING);
    auto magnetic = zero_length_subtractive_gap_magnetic("twoPieceSet");
    auto records = Logger::getInstance().drainCollected();
    Logger::getInstance().disableCollector();
    bool warned = false;
    for (const auto& record : records) {
        if (record.moduleOfOrigin == "Core" && record.message.find("subtractive gap of length zero") != std::string::npos) {
            warned = true;
        }
    }
    CHECK(warned);
    // Every gap is the residual gap of mating halves now.
    for (const auto& gap : magnetic.get_core().get_functional_description().get_gapping()) {
        CHECK(gap.get_type() == GapType::RESIDUAL);
        CHECK(gap.get_length() == Constants().residualGap);
    }
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(100e3, 100e-6, 25, WaveformLabel::TRIANGULAR, 0.4, 0.5, 1.0);
    OpenMagnetics::Mas mas;
    mas.set_inputs(inputs);
    mas.set_magnetic(magnetic);
    MagneticSimulator magneticSimulator;
    OpenMagnetics::Mas simulated;
    REQUIRE_NOTHROW(simulated = magneticSimulator.simulate(mas, true));
    REQUIRE(simulated.get_outputs().size() == 1);
    CHECK(simulated.get_outputs()[0].get_winding_losses());
    settings.reset();
}

TEST_CASE("Test_Zero_Length_Subtractive_Gap_Outside_A_Two_Piece_Set_Is_Refused",
          "[simulator][gapping][zero-length-gap]") {
    settings.reset();
    REQUIRE_THROWS_WITH(zero_length_subtractive_gap_magnetic("pieceAndPlate"), Catch::Matchers::ContainsSubstring("subtractive gap of length zero"));
    settings.reset();
}
