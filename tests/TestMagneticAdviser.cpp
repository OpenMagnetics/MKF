#include "Painter.h"
#include "MagneticAdviser.h"
#include "InputsWrapper.h"
#include "TestingUtils.h"
#include "Settings.h"

#include <UnitTest++.h>
#include <vector>

SUITE(MagneticAdviser) {
    bool plot = false;
    auto settings = OpenMagnetics::Settings::GetInstance();

    TEST(Test_MagneticAdviser_High_Current) {
        srand (time(NULL));

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {24, 78, 76};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 507026;
        double magnetizingInductance = 100e-6;
        double temperature = 25;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::TRIANGULAR;
        double peakToPeak = 100;
        double dutyCycle = 0.5;
        double dcCurrent = 0;

        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        OpenMagnetics::MasWrapper masMagnetic;
        inputs.process_waveforms();

        OpenMagnetics::MagneticAdviser magneticAdviser;
        auto masMagnetics = magneticAdviser.get_advised_magnetic(inputs, 5);

        for (auto [masMagnetic, scoring] : masMagnetics) {
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            if (plot) {
                auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
                auto outFile = outputFilePath;
                std::string filename = "MagneticAdviser" + std::to_string(scoring) + ".svg";
                outFile.append(filename);
                OpenMagnetics::Painter painter(outFile);

                painter.paint_core(masMagnetic.get_mutable_magnetic());
                painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
                painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
                painter.export_svg();
            }
        }
    }

    TEST(Test_MagneticAdviser) {
        srand (time(NULL));

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {24, 78, 76};
        std::vector<int64_t> numberParallels = {1, 1, 1};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 507026;
        double magnetizingInductance = 100e-6;
        double temperature = 25;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::TRIANGULAR;
        double peakToPeak = 1;
        double dutyCycle = 0.5;
        double dcCurrent = 0;

        settings->set_coil_allow_margin_tape(true);
        settings->set_coil_allow_insulated_wire(false);
        settings->set_coil_fill_sections_with_margin_tape(true);

        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);


        auto requirements = inputs.get_mutable_design_requirements();
        auto insulationRequirements = requirements.get_insulation().value();
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
        insulationRequirements.set_standards(standards);
        requirements.set_insulation(insulationRequirements);
        inputs.set_design_requirements(requirements);

        OpenMagnetics::MasWrapper masMagnetic;
        inputs.process_waveforms();

        OpenMagnetics::MagneticAdviser magneticAdviser;
        auto masMagnetics = magneticAdviser.get_advised_magnetic(inputs, 4);
        CHECK(masMagnetics.size() > 0);

        auto scorings = magneticAdviser.get_scorings();
        // std::cout << "scorings.size(): " << scorings.size() << std::endl;
        for (auto [name, values] : scorings) {
            double scoringTotal = 0;
            // std::cout << name << std::endl;
            for (auto [key, value] : values) {
                // std::cout << magic_enum::enum_name(key) << std::endl;
                // std::cout << value << std::endl;
                scoringTotal += value;
            }
            // std::cout << "scoringTotal: " << scoringTotal << std::endl;
        }

        for (auto masMagneticWithScoring : masMagnetics) {
            auto masMagnetic = masMagneticWithScoring.first;
            // std::cout << "name: " << masMagnetic.get_magnetic().get_manufacturer_info().value().get_reference().value() << std::endl;
            // std::cout << "scoringSecond: " << masMagneticWithScoring.second << std::endl;
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            if (plot) {
                auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
                auto outFile = outputFilePath;
                std::string filename = "MagneticAdviser" + masMagnetic.get_magnetic().get_manufacturer_info().value().get_reference().value() + ".svg";
                outFile.append(filename);
                OpenMagnetics::Painter painter(outFile);

                painter.paint_core(masMagnetic.get_mutable_magnetic());
                painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
                painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
                painter.export_svg();
            }
        }
    }

    TEST(Test_MagneticAdviser_No_Insulation_Requirements) {
        srand (time(NULL));

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {24, 78, 76};
        std::vector<int64_t> numberParallels = {1, 1, 1};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 507026;
        double magnetizingInductance = 100e-6;
        double temperature = 25;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::TRIANGULAR;
        double peakToPeak = 1;
        double dutyCycle = 0.5;
        double dcCurrent = 0;

        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->set_coil_allow_margin_tape(true);
        settings->set_coil_allow_insulated_wire(false);
        settings->set_coil_fill_sections_with_margin_tape(true);

        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);


        auto requirements = inputs.get_mutable_design_requirements();
        requirements.set_insulation(std::nullopt);
        inputs.set_design_requirements(requirements);

        OpenMagnetics::MasWrapper masMagnetic;
        inputs.process_waveforms();

        OpenMagnetics::MagneticAdviser magneticAdviser;
        auto masMagnetics = magneticAdviser.get_advised_magnetic(inputs, 1);
        CHECK(masMagnetics.size() > 0);

        for (auto masMagneticWithScoring : masMagnetics) {
            auto masMagnetic = masMagneticWithScoring.first;
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            if (plot) {
                auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
                auto outFile = outputFilePath;
                std::string filename = "MagneticAdviser" + std::to_string(std::rand()) + ".svg";
                outFile.append(filename);
                OpenMagnetics::Painter painter(outFile);

                painter.paint_core(masMagnetic.get_mutable_magnetic());
                painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
                painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
                painter.export_svg();
            }
        }

    }

    TEST(MagneticAdviserJsonHV) {
        auto settings = OpenMagnetics::Settings::GetInstance();
        OpenMagnetics::clear_databases();
        settings->reset();
        settings->set_use_only_cores_in_stock(true);
        settings->set_use_toroidal_cores(true);
        settings->set_use_concentric_cores(true);
        settings->set_coil_fill_sections_with_margin_tape(true);
        std::string masString = R"({"inputs": {"designRequirements": {"insulation": {"altitude": {"maximum": 1000 }, "cti": "Group IIIB", "pollutionDegree": "P2", "overvoltageCategory": "OVC-II", "insulationType": "Basic", "mainSupplyVoltage": {"nominal": null, "minimum": 100, "maximum": 815 }, "standards": ["IEC 60664-1"] }, "leakageInductance": [{"maximum": 0.0000027 }, {"maximum": 7.2e-7 } ], "magnetizingInductance": {"maximum": null, "minimum": null, "nominal": 0.0007 }, "market": "Commercial", "maximumDimensions": {"width": null, "height": 0.02, "depth": null }, "maximumWeight": null, "name": "My Design Requirements", "operatingTemperature": null, "strayCapacitance": [{"maximum": 5e-11 }, {"maximum": 5e-11 } ], "terminalType": ["Pin", "Pin", "Pin"], "topology": "Flyback Converter", "turnsRatios": [{"nominal": 3.6 }, {"nominal": 7 } ] }, "operatingPoints": [{"name": "Operating Point No. 1", "conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"name": "Primary winding excitation", "frequency": 42000, "current": {"waveform": {"ancillaryLabel": null, "data": [0, 0, 1.1, 0, 0 ], "numberPeriods": null, "time": [0, 0, 0.000007, 0.000007, 0.00002380952380952381 ] }, "processed": {"dutyCycle": 0.294, "peakToPeak": 1.1, "offset": 0, "label": "Flyback Primary", "acEffectiveFrequency": 320324.70197566116, "effectiveFrequency": 299906.6380380921, "peak": 1.0815263605442182, "rms": 0.34251379452955516, "thd": 1.070755363168311 }, "harmonics": {"amplitudes": [0.16053906914328236, 0.29202196202762337, 0.21752279286410783, 0.1310845319633303, 0.07713860380309742, 0.06991594662329052, 0.06385479128189217, 0.04886550825820873, 0.041656942163512685, 0.04106933856014564, 0.036071164493133245, 0.03055299653431604, 0.02975135660550652, 0.028391215457215468, 0.02482068925428346, 0.023301839916889198, 0.023120564125526488, 0.021217296242203084, 0.019422484452150467, 0.019302597014669035, 0.018564086064230358, 0.016990298638197066, 0.016530128240448857, 0.016408746965119517, 0.015342013775419296, 0.014571250724490458, 0.014600261585603161, 0.014080719038131833, 0.013228104677648548, 0.013112646237623088, 0.012998033401152062, 0.012300975273665412, 0.011953654437324483, 0.012016309039770377, 0.011608286047010533, 0.011111174591424324, 0.011137874398755002, 0.01101797292674616, 0.010530630511366074, 0.01040049906973778, 0.01046208177619542, 0.010124703775041825, 0.009838808475080322, 0.009931053977331523, 0.0098034018119726, 0.009456464631666566, 0.009454698580407995, 0.009502890216301076, 0.009218272478811467, 0.00907613407579727, 0.009199516076645192, 0.00906471876755199, 0.008824889331959554, 0.008907818118804589, 0.008937898630441313, 0.00869876368067147, 0.008666288051213203, 0.008803433740422266, 0.008663604885115571, 0.008515486625290338, 0.008660125560530002, 0.008670306514975606, 0.008475927879718145, 0.008536460912837853 ], "frequencies": [0, 42000, 84000, 126000, 168000, 210000, 252000, 294000, 336000, 378000, 420000, 462000, 504000, 546000, 588000, 630000, 672000, 714000, 756000, 798000, 840000, 882000, 924000, 966000, 1008000, 1050000, 1092000, 1134000, 1176000, 1218000, 1260000, 1302000, 1344000, 1386000, 1428000, 1470000, 1512000, 1554000, 1596000, 1638000, 1680000, 1722000, 1764000, 1806000, 1848000, 1890000, 1932000, 1974000, 2016000, 2058000, 2100000, 2142000, 2184000, 2226000, 2268000, 2310000, 2352000, 2394000, 2436000, 2478000, 2520000, 2562000, 2604000, 2646000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-29.4, 70.6, 70.6, -29.4, -29.4 ], "numberPeriods": null, "time": [0, 0, 0.000007, 0.000007, 0.00002380952380952381 ] }, "processed": {"dutyCycle": 0.294, "peakToPeak": 100, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 273993.1760953399, "effectiveFrequency": 273985.04815560044, "peak": 70.6, "rms": 45.33538904652732, "thd": 0.7943280084778311 }, "harmonics": {"amplitudes": [0.4937500000000009, 50.19273146886441, 30.88945734421223, 8.60726309175931, 7.514577540357868, 12.576487416046705, 7.890222788854489, 0.6723406687769768, 6.659332230305735, 6.771471403422372, 2.1663923211140617, 3.134209960512532, 5.356728264586365, 3.4347141143172375, 0.6805380572551197, 3.777520625726094, 3.7722086912079633, 1.0283716619016732, 2.1769829843489013, 3.474177141872247, 2.1027712520453545, 0.6945565902694819, 2.7474727308098412, 2.608727438661375, 0.5486770349329464, 1.7724349572072589, 2.6198081625264167, 1.4623766283543616, 0.71496628279737, 2.2318766983434566, 1.995656440067472, 0.2773827898281461, 1.562500000000001, 2.1411881515502134, 1.0841275629577949, 0.7426431917434089, 1.9342807953981735, 1.620450382955772, 0.09545259242683589, 1.4471144811333827, 1.8430944173488684, 0.8315259963159228, 0.7788650462348548, 1.7521633971823452, 1.3695427043271637, 0.042929496561759636, 1.388303796831737, 1.6472150875774871, 0.6472086912079588, 0.8254655623504514, 1.6415464278026473, 1.1917877639027972, 0.16004318769545572, 1.3696996040394052, 1.5166136111631239, 0.5023333497836172, 0.8850851521101957, 1.581550490460967, 1.0607923422947243, 0.2691521002455623, 1.3846695378376832, 1.4323337135341552, 0.38011439553529425, 0.9615889729028446 ], "frequencies": [0, 42000, 84000, 126000, 168000, 210000, 252000, 294000, 336000, 378000, 420000, 462000, 504000, 546000, 588000, 630000, 672000, 714000, 756000, 798000, 840000, 882000, 924000, 966000, 1008000, 1050000, 1092000, 1134000, 1176000, 1218000, 1260000, 1302000, 1344000, 1386000, 1428000, 1470000, 1512000, 1554000, 1596000, 1638000, 1680000, 1722000, 1764000, 1806000, 1848000, 1890000, 1932000, 1974000, 2016000, 2058000, 2100000, 2142000, 2184000, 2226000, 2268000, 2310000, 2352000, 2394000, 2436000, 2478000, 2520000, 2562000, 2604000, 2646000 ] } } }, {"name": "Primary winding excitation", "frequency": 42000, "current": {"waveform": {"ancillaryLabel": null, "data": [0, 0, 3.96, 0, 0 ], "numberPeriods": null, "time": [0, 0.000011904761904761905, 0.000011904761904761905, 0.00002380952380952381, 0.00002380952380952381 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 3.96, "offset": 0, "label": "Flyback Secondary", "acEffectiveFrequency": 273605.22552591236, "effectiveFrequency": 239619.35667760254, "peak": 3.9599999999999986, "rms": 1.635596311125932, "thd": 0.6765065136547783 }, "harmonics": {"amplitudes": [1.0054687500000032, 1.5109819803933975, 0.6305067526445318, 0.43631079159665215, 0.3156335707813906, 0.25867345007829645, 0.2108457708927511, 0.1845098972167098, 0.15858039332900484, 0.14374429235357594, 0.12732519999086678, 0.11798991339580255, 0.1065764142034965, 0.10027664365108102, 0.09183268563443873, 0.08737513814977015, 0.08084358345172529, 0.0775829517661027, 0.07235909588115942, 0.06991721285428486, 0.06562943182065245, 0.06377074179801741, 0.060177659717708575, 0.05874824283059186, 0.05568602880656547, 0.054581360149729106, 0.051934759310085545, 0.0510815825974439, 0.0487670673260459, 0.048112679051860716, 0.04606818176826717, 0.04557382535330518, 0.043752232085917345, 0.04338887874386873, 0.04175376612752825, 0.04149933709574444, 0.04002209784539918, 0.03985958700992044, 0.038517442612831645, 0.03843361846308103, 0.03720821487910577, 0.037192705080064635, 0.03606910115775439, 0.036113735798427, 0.03507966215426201, 0.03517799562893598, 0.03422330356713118, 0.0343702622545533, 0.03348650869654576, 0.03367812888504953, 0.03285826032368028, 0.03309149205770409, 0.032329601799393176, 0.032602161771072315, 0.0318933023014623, 0.032203563952980456, 0.031543601457069734, 0.03189051394290526, 0.03127601564366044, 0.03165904576198601, 0.03108719333289207, 0.0315062863179612, 0.03097481051580278, 0.03143036691725387 ], "frequencies": [0, 42000, 84000, 126000, 168000, 210000, 252000, 294000, 336000, 378000, 420000, 462000, 504000, 546000, 588000, 630000, 672000, 714000, 756000, 798000, 840000, 882000, 924000, 966000, 1008000, 1050000, 1092000, 1134000, 1176000, 1218000, 1260000, 1302000, 1344000, 1386000, 1428000, 1470000, 1512000, 1554000, 1596000, 1638000, 1680000, 1722000, 1764000, 1806000, 1848000, 1890000, 1932000, 1974000, 2016000, 2058000, 2100000, 2142000, 2184000, 2226000, 2268000, 2310000, 2352000, 2394000, 2436000, 2478000, 2520000, 2562000, 2604000, 2646000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-13.89, 13.89, 13.89, -13.89, -13.89 ], "numberPeriods": null, "time": [0, 0, 0.000011904761904761905, 0.000011904761904761905, 0.00002380952380952381 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 27.78, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 248423.92512497256, "effectiveFrequency": 248408.7565140517, "peak": 13.89, "rms": 13.890000000000022, "thd": 0.48331514845248535 }, "harmonics": {"amplitudes": [0.21703125, 17.681745968226167, 0.4340625, 5.884441743008607, 0.4340625, 3.5192857754085134, 0.4340625, 2.50156382659689, 0.4340625, 1.932968090534883, 0.4340625, 1.56850033166751, 0.4340625, 1.313925940874188, 0.4340625, 1.1252647178556892, 0.4340625, 0.9792293094780014, 0.4340625, 0.8623340820515439, 0.4340625, 0.7662274695502621, 0.4340625, 0.6854595927802318, 0.4340625, 0.6163213952979067, 0.4340625, 0.5561996920846202, 0.4340625, 0.5031990666519233, 0.4340625, 0.45591010107099894, 0.4340625, 0.41326185461487136, 0.4340625, 0.3744248874701925, 0.4340625, 0.33874569976854113, 0.4340625, 0.30570130348173624, 0.4340625, 0.2748670467095772, 0.4340625, 0.24589336899764186, 0.4340625, 0.21848870156913297, 0.4340625, 0.1924066733732579, 0.4340625, 0.16743638267205946, 0.4340625, 0.1433948809785266, 0.4340625, 0.12012127132032469, 0.4340625, 0.09747199388796957, 0.4340625, 0.07531698847858959, 0.4340625, 0.05353650312310587, 0.4340625, 0.03201837355771753, 0.4340625, 0.0106556362841701 ], "frequencies": [0, 42000, 84000, 126000, 168000, 210000, 252000, 294000, 336000, 378000, 420000, 462000, 504000, 546000, 588000, 630000, 672000, 714000, 756000, 798000, 840000, 882000, 924000, 966000, 1008000, 1050000, 1092000, 1134000, 1176000, 1218000, 1260000, 1302000, 1344000, 1386000, 1428000, 1470000, 1512000, 1554000, 1596000, 1638000, 1680000, 1722000, 1764000, 1806000, 1848000, 1890000, 1932000, 1974000, 2016000, 2058000, 2100000, 2142000, 2184000, 2226000, 2268000, 2310000, 2352000, 2394000, 2436000, 2478000, 2520000, 2562000, 2604000, 2646000 ] } } }, {"name": "Primary winding excitation", "frequency": 42000, "current": {"waveform": {"ancillaryLabel": null, "data": [0, 0, 0.08, 0, 0 ], "numberPeriods": null, "time": [0, 0.000011904761904761905, 0.000011904761904761905, 0.00002380952380952381, 0.00002380952380952381 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 0.08, "offset": 0, "label": "Flyback Secondary", "acEffectiveFrequency": 273605.2255259123, "effectiveFrequency": 239619.35667760245, "peak": 0.07999999999999997, "rms": 0.03304234971971578, "thd": 0.6765065136547783 }, "harmonics": {"amplitudes": [0.020312500000000067, 0.030524888492795912, 0.012737510154434986, 0.008814359426194993, 0.006376435773361427, 0.00522572626420801, 0.004259510523085881, 0.0037274726710446416, 0.0032036443096768643, 0.0029039250980520386, 0.0025722262624417523, 0.002383634614056616, 0.002153058872797908, 0.0020257907808299194, 0.0018552057703927027, 0.0017651543060559615, 0.0016332037060954627, 0.0015673323589111663, 0.0014617999167910994, 0.0014124689465512088, 0.001325847107487929, 0.0012882978141013626, 0.0012157102973274464, 0.001186833188496806, 0.001124970278920515, 0.0011026537403985668, 0.0010491870567694052, 0.0010319511635847236, 0.0009851932793140594, 0.000971973314179004, 0.0009306703387528712, 0.0009206833404708095, 0.0008838834764831795, 0.0008765430049266404, 0.0008435104268187535, 0.0008383704463786768, 0.000808527229199984, 0.0008052441820185946, 0.0007781301537945788, 0.0007764367366278998, 0.0007516811086688038, 0.0007513677793952459, 0.0007286687102576646, 0.0007295704201702419, 0.0007086800435204446, 0.0007106665783623426, 0.000691379870043054, 0.0006943487324152177, 0.0006764951251827419, 0.00068036624010201, 0.0006638032388622281, 0.0006685149910647295, 0.000653123268674609, 0.0006586295307287333, 0.0006443091374032789, 0.0006505770495551611, 0.0006372444738801973, 0.0006442528069273796, 0.0006318386998719282, 0.0006395766820603246, 0.0006280241077351936, 0.0006364906326860852, 0.0006257537477939952, 0.0006349569074192705 ], "frequencies": [0, 42000, 84000, 126000, 168000, 210000, 252000, 294000, 336000, 378000, 420000, 462000, 504000, 546000, 588000, 630000, 672000, 714000, 756000, 798000, 840000, 882000, 924000, 966000, 1008000, 1050000, 1092000, 1134000, 1176000, 1218000, 1260000, 1302000, 1344000, 1386000, 1428000, 1470000, 1512000, 1554000, 1596000, 1638000, 1680000, 1722000, 1764000, 1806000, 1848000, 1890000, 1932000, 1974000, 2016000, 2058000, 2100000, 2142000, 2184000, 2226000, 2268000, 2310000, 2352000, 2394000, 2436000, 2478000, 2520000, 2562000, 2604000, 2646000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-7.14, 7.14, 7.14, -7.14, -7.14 ], "numberPeriods": null, "time": [0, 0, 0.000011904761904761905, 0.000011904761904761905, 0.00002380952380952381 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 14.28, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 248423.9251249721, "effectiveFrequency": 248408.75651405123, "peak": 7.14, "rms": 7.139999999999993, "thd": 0.4833151484524845 }, "harmonics": {"amplitudes": [0.1115625, 9.08910483895859, 0.223125, 3.0248318246998886, 0.223125, 1.8090497074454124, 0.223125, 1.2859010598921377, 0.223125, 0.9936207463224668, 0.223125, 0.8062701488917218, 0.223125, 0.6754090149634054, 0.223125, 0.5784298117703113, 0.223125, 0.5033619344616945, 0.223125, 0.44327324304161464, 0.223125, 0.3938707078897675, 0.223125, 0.3523528792261233, 0.223125, 0.3168131578421202, 0.223125, 0.285908265045658, 0.223125, 0.2586638830737746, 0.223125, 0.23435551631727392, 0.223125, 0.21243265960764457, 0.223125, 0.19246894863478564, 0.223125, 0.1741284590602866, 0.223125, 0.1571423547055144, 0.223125, 0.14129234798462062, 0.223125, 0.1263987512342087, 0.223125, 0.1123116867677183, 0.223125, 0.09890451028690123, 0.223125, 0.08606881009924447, 0.223125, 0.07371054356995504, 0.223125, 0.0617470034000806, 0.223125, 0.05010439426638574, 0.223125, 0.0387158601682599, 0.223125, 0.027519843938011213, 0.223125, 0.016458688783449027, 0.223125, 0.005477411308063118 ], "frequencies": [0, 42000, 84000, 126000, 168000, 210000, 252000, 294000, 336000, 378000, 420000, 462000, 504000, 546000, 588000, 630000, 672000, 714000, 756000, 798000, 840000, 882000, 924000, 966000, 1008000, 1050000, 1092000, 1134000, 1176000, 1218000, 1260000, 1302000, 1344000, 1386000, 1428000, 1470000, 1512000, 1554000, 1596000, 1638000, 1680000, 1722000, 1764000, 1806000, 1848000, 1890000, 1932000, 1974000, 2016000, 2058000, 2100000, 2142000, 2184000, 2226000, 2268000, 2310000, 2352000, 2394000, 2436000, 2478000, 2520000, 2562000, 2604000, 2646000 ] } } } ] } ] }, "magnetic": {"coil": {"functionalDescription": [{"name": "Primary"}, {"name": "Secondary"}, {"name": "Tertiary"} ] } } })";
        json masJson = json::parse(masString);

        if (masJson["inputs"]["designRequirements"]["insulation"]["cti"] == "GROUP_IIIa") {
            masJson["inputs"]["designRequirements"]["insulation"]["cti"] = "GROUP_IIIA";
        }
        if (masJson["inputs"]["designRequirements"]["insulation"]["cti"] == "GROUP_IIIb") {
            masJson["inputs"]["designRequirements"]["insulation"]["cti"] = "GROUP_IIIB";
        }

        OpenMagnetics::InputsWrapper inputs(masJson["inputs"]);
        OpenMagnetics::MasWrapper masMagnetic;

        OpenMagnetics::MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 2);
        CHECK(masMagnetics.size() > 0);

        for (auto [masMagnetic, scoring] : masMagnetics) {
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "MagneticAdviser_MagneticAdviserJsonHV_" + std::to_string(scoring) + ".svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagnetic.get_mutable_magnetic());
            painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            painter.export_svg();
        }
    }

    TEST(MagneticAdviserJsonLV) {
        auto settings = OpenMagnetics::Settings::GetInstance();
        OpenMagnetics::clear_databases();
        settings->reset();
        settings->set_use_only_cores_in_stock(true);
        settings->set_coil_fill_sections_with_margin_tape(true);
        std::string masString = R"({"inputs": {"designRequirements": {"insulation": {"altitude": {"maximum": 1000 }, "cti": "Group IIIB", "pollutionDegree": "P2", "overvoltageCategory": "OVC-II", "insulationType": "Basic", "mainSupplyVoltage": {"nominal": 120, "minimum": 20, "maximum": 450 }, "standards": ["IEC 60664-1"] }, "leakageInductance": [{"maximum": 0.00000135 }, {"maximum": 0.00000135 } ], "magnetizingInductance": {"maximum": 0.000088, "minimum": 0.000036, "nominal": 0.000039999999999999996 }, "market": "Commercial", "maximumDimensions": {"width": null, "height": 0.02, "depth": null }, "maximumWeight": null, "name": "My Design Requirements", "operatingTemperature": null, "strayCapacitance": [{"maximum": 5e-11 }, {"maximum": 5e-11 } ], "terminalType": ["Pin", "Pin", "Pin"], "topology": "Flyback Converter", "turnsRatios": [{"nominal": 1.2 }, {"nominal": 1.2 } ] }, "operatingPoints": [{"name": "Operating Point No. 1", "conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"name": "Primary winding excitation", "frequency": 36000, "current": {"waveform": {"ancillaryLabel": null, "data": [0, 0, 5.5, 0, 0 ], "numberPeriods": null, "time": [0, 0, 0.000011111111111111112, 0.000011111111111111112, 0.00002777777777777778 ] }, "processed": {"dutyCycle": 0.4, "peakToPeak": 5.5, "offset": 0, "label": "Flyback Primary", "acEffectiveFrequency": 247671.20334838802, "effectiveFrequency": 224591.09893738205, "peak": 5.478515625000006, "rms": 2.025897603604934, "thd": 0.8123116860415294 }, "harmonics": {"amplitudes": [1.1128234863281259, 1.8581203560604207, 1.051281221517565, 0.5159082238421047, 0.4749772671818343, 0.35101798151542746, 0.28680520988943464, 0.26606329358167974, 0.2109058168232097, 0.2038324242201338, 0.17683391011628416, 0.15977271740890614, 0.15330090420784062, 0.13359870506402888, 0.1311678604251054, 0.11938227872282572, 0.11196939263234992, 0.10881687126769939, 0.09881440376099933, 0.09781517962337143, 0.09113803697574208, 0.08720137281111257, 0.08534603089826898, 0.07933100843010468, 0.07894777491614259, 0.0746052794656816, 0.07229734773884018, 0.07109041486606842, 0.0670974294019412, 0.06703716258861198, 0.06396085871081211, 0.0625474905371643, 0.06171647321147238, 0.05889230550248535, 0.059028047581821314, 0.05671822612277054, 0.05585484049699437, 0.055265364667603324, 0.0531817689812044, 0.053452567468547284, 0.05164394009764037, 0.051150939884780554, 0.05072953903613357, 0.049149522805038266, 0.049525208170729855, 0.04806384343218855, 0.04784167797419847, 0.047545488047249944, 0.04632889186784806, 0.046795775944523485, 0.045586400651502844, 0.045578916197115495, 0.045382141359988915, 0.04444264358178238, 0.04499725549277023, 0.04397809673849301, 0.04415534166549457, 0.044042673430388084, 0.04332652563384539, 0.04397294643650677, 0.043102730085842, 0.043452676225114936, 0.04341597260904029, 0.04289100286968138 ], "frequencies": [0, 36000, 72000, 108000, 144000, 180000, 216000, 252000, 288000, 324000, 360000, 396000, 432000, 468000, 504000, 540000, 576000, 612000, 648000, 684000, 720000, 756000, 792000, 828000, 864000, 900000, 936000, 972000, 1008000, 1044000, 1080000, 1116000, 1152000, 1188000, 1224000, 1260000, 1296000, 1332000, 1368000, 1404000, 1440000, 1476000, 1512000, 1548000, 1584000, 1620000, 1656000, 1692000, 1728000, 1764000, 1800000, 1836000, 1872000, 1908000, 1944000, 1980000, 2016000, 2052000, 2088000, 2124000, 2160000, 2196000, 2232000, 2268000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-8, 12, 12, -8, -8 ], "numberPeriods": null, "time": [0, 0, 0.000011111111111111112, 0.000011111111111111112, 0.00002777777777777778 ] }, "processed": {"dutyCycle": 0.4, "peakToPeak": 20, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 217445.61780293446, "effectiveFrequency": 217445.06394458748, "peak": 12, "rms": 9.79157801378307, "thd": 0.5579294461931013 }, "harmonics": {"amplitudes": [0.03125, 12.090982167347164, 3.7938629699811433, 2.4460154592542187, 3.050934294867692, 0.06265085868474476, 2.0052584267980813, 1.1245773748206938, 0.8899247078195547, 1.3746157728063957, 0.06310658027096225, 1.0931619934351904, 0.682943044655975, 0.5329821111108219, 0.8998037891938968, 0.06387675236155328, 0.7544417382415929, 0.5036912732460003, 0.37575767434186663, 0.6781566555968945, 0.06497787099509279, 0.579613707579126, 0.4082105683349906, 0.28790910341678466, 0.5516771452281819, 0.06643416431010025, 0.4742269329016237, 0.35024440150573616, 0.23220844850577863, 0.4713901788294872, 0.0682788501833682, 0.40480954492681936, 0.3125, 0.19399945657286152, 0.41719035005974064, 0.07055595095604376, 0.35652867811806305, 0.28709106236926035, 0.16634675530328177, 0.3793488450042862, 0.07332285477181877, 0.32183350571658953, 0.2699539504381825, 0.14554626367323475, 0.3526337777044684, 0.0766539126369276, 0.29650814745710025, 0.25883312123614344, 0.12944173824159266, 0.33403609127672834, 0.08064551519955634, 0.27805045130481554, 0.25243555700995646, 0.1166943265232067, 0.3217665194149353, 0.08542334196565987, 0.2649247078195551, 0.250039428516803, 0.10643002348369442, 0.31477312120987916, 0.09115288771096393, 0.25618941041585475, 0.2513050615401854, 0.098055076230366 ], "frequencies": [0, 36000, 72000, 108000, 144000, 180000, 216000, 252000, 288000, 324000, 360000, 396000, 432000, 468000, 504000, 540000, 576000, 612000, 648000, 684000, 720000, 756000, 792000, 828000, 864000, 900000, 936000, 972000, 1008000, 1044000, 1080000, 1116000, 1152000, 1188000, 1224000, 1260000, 1296000, 1332000, 1368000, 1404000, 1440000, 1476000, 1512000, 1548000, 1584000, 1620000, 1656000, 1692000, 1728000, 1764000, 1800000, 1836000, 1872000, 1908000, 1944000, 1980000, 2016000, 2052000, 2088000, 2124000, 2160000, 2196000, 2232000, 2268000 ] } } }, {"name": "Primary winding excitation", "frequency": 36000, "current": {"waveform": {"ancillaryLabel": null, "data": [0, 0, 6.6, 0, 0 ], "numberPeriods": null, "time": [0, 0.00001388888888888889, 0.00001388888888888889, 0.00002777777777777778, 0.00002777777777777778 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 6.6, "offset": 0, "label": "Flyback Secondary", "acEffectiveFrequency": 234518.76473649617, "effectiveFrequency": 205388.02000937346, "peak": 6.599999999999992, "rms": 2.7259938518765505, "thd": 0.6765065136547781 }, "harmonics": {"amplitudes": [1.6757812500000036, 2.5183033006556617, 1.0508445877408874, 0.7271846526610849, 0.5260559513023173, 0.43112241679715924, 0.35140961815458516, 0.30751649536118214, 0.2643006555483416, 0.23957382058929272, 0.21220866665144478, 0.19664985565967041, 0.17762735700582763, 0.16712773941846806, 0.15305447605739808, 0.14562523024961643, 0.13473930575287565, 0.12930491961017093, 0.12059849313526576, 0.11652868809047452, 0.10938238636775419, 0.10628456966336206, 0.10029609952951436, 0.0979137380509861, 0.09281004801094254, 0.09096893358288156, 0.08655793218347599, 0.08513597099573951, 0.08127844554340995, 0.08018779841976777, 0.07678030294711186, 0.07595637558884157, 0.07292038680986246, 0.07231479790644763, 0.06958961021254718, 0.06916556182624073, 0.06670349640899874, 0.06643264501653386, 0.06419573768805274, 0.06405603077180153, 0.062013691465176275, 0.0619878418001075, 0.060115168596257366, 0.060189559664044844, 0.058466103590436726, 0.05862999271489314, 0.05703883927855185, 0.05728377042425536, 0.05581084782757627, 0.05613021480841572, 0.054763767206133816, 0.05515248676284005, 0.05388266966565528, 0.05433693628512042, 0.05315550383577058, 0.05367260658830061, 0.05257266909511635, 0.05315085657150875, 0.052126692739434176, 0.05276507626997662, 0.05181198888815346, 0.05251047719660189, 0.051624684193004804, 0.05238394486208969 ], "frequencies": [0, 36000, 72000, 108000, 144000, 180000, 216000, 252000, 288000, 324000, 360000, 396000, 432000, 468000, 504000, 540000, 576000, 612000, 648000, 684000, 720000, 756000, 792000, 828000, 864000, 900000, 936000, 972000, 1008000, 1044000, 1080000, 1116000, 1152000, 1188000, 1224000, 1260000, 1296000, 1332000, 1368000, 1404000, 1440000, 1476000, 1512000, 1548000, 1584000, 1620000, 1656000, 1692000, 1728000, 1764000, 1800000, 1836000, 1872000, 1908000, 1944000, 1980000, 2016000, 2052000, 2088000, 2124000, 2160000, 2196000, 2232000, 2268000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-8.335, 8.335, 8.335, -8.335, -8.335 ], "numberPeriods": null, "time": [0, 0, 0.00001388888888888889, 0.00001388888888888889, 0.00002777777777777778 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 16.67, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 212934.79296426187, "effectiveFrequency": 212921.79129775826, "peak": 8.335, "rms": 8.33500000000001, "thd": 0.48331514845248513 }, "harmonics": {"amplitudes": [0.130234375, 10.610320564806702, 0.26046875, 3.531088691718988, 0.26046875, 2.1118248335514735, 0.26046875, 1.5011183941457937, 0.26046875, 1.1599200168904427, 0.26046875, 0.9412131219905466, 0.26046875, 0.7884501596246478, 0.26046875, 0.6752398432920929, 0.26046875, 0.587608084557174, 0.26046875, 0.5174625323181874, 0.26046875, 0.4597916456948479, 0.26046875, 0.4113251048108877, 0.26046875, 0.36983720876947807, 0.26046875, 0.333759858425148, 0.26046875, 0.3019556674257579, 0.26046875, 0.2735788835440449, 0.26046875, 0.24798686524225744, 0.26046875, 0.2246818889174983, 0.26046875, 0.20327180760048874, 0.26046875, 0.18344279082219334, 0.26046875, 0.16494001687000165, 0.26046875, 0.14755372430492025, 0.26046875, 0.1311089508695984, 0.26046875, 0.1154578561962637, 0.26046875, 0.10047388405843227, 0.26046875, 0.0860472521926578, 0.26046875, 0.0720814108318868, 0.26046875, 0.058490213754947895, 0.26046875, 0.04519561547653339, 0.26046875, 0.032125756193742694, 0.26046875, 0.019213329273115587, 0.26046875, 0.006394148914942832 ], "frequencies": [0, 36000, 72000, 108000, 144000, 180000, 216000, 252000, 288000, 324000, 360000, 396000, 432000, 468000, 504000, 540000, 576000, 612000, 648000, 684000, 720000, 756000, 792000, 828000, 864000, 900000, 936000, 972000, 1008000, 1044000, 1080000, 1116000, 1152000, 1188000, 1224000, 1260000, 1296000, 1332000, 1368000, 1404000, 1440000, 1476000, 1512000, 1548000, 1584000, 1620000, 1656000, 1692000, 1728000, 1764000, 1800000, 1836000, 1872000, 1908000, 1944000, 1980000, 2016000, 2052000, 2088000, 2124000, 2160000, 2196000, 2232000, 2268000 ] } } }, {"name": "Primary winding excitation", "frequency": 36000, "current": {"waveform": {"ancillaryLabel": null, "data": [0, 0, 0.08, 0, 0 ], "numberPeriods": null, "time": [0, 0.00001388888888888889, 0.00001388888888888889, 0.00002777777777777778, 0.00002777777777777778 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 0.08, "offset": 0, "label": "Flyback Secondary", "acEffectiveFrequency": 234518.7647364962, "effectiveFrequency": 205388.02000937346, "peak": 0.0799999999999999, "rms": 0.033042349719715765, "thd": 0.6765065136547783 }, "harmonics": {"amplitudes": [0.020312500000000046, 0.030524888492795898, 0.012737510154435002, 0.008814359426194969, 0.006376435773361423, 0.005225726264207993, 0.004259510523085882, 0.003727472671044633, 0.0032036443096768678, 0.0029039250980520334, 0.002572226262441755, 0.002383634614056611, 0.0021530588727979106, 0.0020257907808299163, 0.001855205770392704, 0.0017651543060559591, 0.0016332037060954637, 0.0015673323589111632, 0.0014617999167910998, 0.0014124689465512055, 0.0013258471074879294, 0.0012882978141013591, 0.0012157102973274472, 0.0011868331884968021, 0.0011249702789205155, 0.0011026537403985642, 0.0010491870567694056, 0.0010319511635847216, 0.0009851932793140603, 0.0009719733141790018, 0.0009306703387528721, 0.0009206833404708108, 0.0008838834764831795, 0.0008765430049266376, 0.0008435104268187531, 0.0008383704463786751, 0.0008085272291999843, 0.0008052441820185928, 0.0007781301537945793, 0.0007764367366278984, 0.0007516811086688035, 0.0007513677793952425, 0.0007286687102576648, 0.0007295704201702394, 0.0007086800435204447, 0.0007106665783623409, 0.0006913798700430537, 0.0006943487324152148, 0.0006764951251827443, 0.0006803662401020098, 0.0006638032388622283, 0.0006685149910647278, 0.0006531232686746098, 0.0006586295307287323, 0.0006443091374032788, 0.0006505770495551569, 0.0006372444738801973, 0.0006442528069273779, 0.0006318386998719287, 0.0006395766820603229, 0.0006280241077351938, 0.0006364906326860843, 0.0006257537477939962, 0.0006349569074192687 ], "frequencies": [0, 36000, 72000, 108000, 144000, 180000, 216000, 252000, 288000, 324000, 360000, 396000, 432000, 468000, 504000, 540000, 576000, 612000, 648000, 684000, 720000, 756000, 792000, 828000, 864000, 900000, 936000, 972000, 1008000, 1044000, 1080000, 1116000, 1152000, 1188000, 1224000, 1260000, 1296000, 1332000, 1368000, 1404000, 1440000, 1476000, 1512000, 1548000, 1584000, 1620000, 1656000, 1692000, 1728000, 1764000, 1800000, 1836000, 1872000, 1908000, 1944000, 1980000, 2016000, 2052000, 2088000, 2124000, 2160000, 2196000, 2232000, 2268000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-8.335, 8.335, 8.335, -8.335, -8.335 ], "numberPeriods": null, "time": [0, 0, 0.00001388888888888889, 0.00001388888888888889, 0.00002777777777777778 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 16.67, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 212934.79296426187, "effectiveFrequency": 212921.79129775826, "peak": 8.335, "rms": 8.33500000000001, "thd": 0.48331514845248513 }, "harmonics": {"amplitudes": [0.130234375, 10.610320564806702, 0.26046875, 3.531088691718988, 0.26046875, 2.1118248335514735, 0.26046875, 1.5011183941457937, 0.26046875, 1.1599200168904427, 0.26046875, 0.9412131219905466, 0.26046875, 0.7884501596246478, 0.26046875, 0.6752398432920929, 0.26046875, 0.587608084557174, 0.26046875, 0.5174625323181874, 0.26046875, 0.4597916456948479, 0.26046875, 0.4113251048108877, 0.26046875, 0.36983720876947807, 0.26046875, 0.333759858425148, 0.26046875, 0.3019556674257579, 0.26046875, 0.2735788835440449, 0.26046875, 0.24798686524225744, 0.26046875, 0.2246818889174983, 0.26046875, 0.20327180760048874, 0.26046875, 0.18344279082219334, 0.26046875, 0.16494001687000165, 0.26046875, 0.14755372430492025, 0.26046875, 0.1311089508695984, 0.26046875, 0.1154578561962637, 0.26046875, 0.10047388405843227, 0.26046875, 0.0860472521926578, 0.26046875, 0.0720814108318868, 0.26046875, 0.058490213754947895, 0.26046875, 0.04519561547653339, 0.26046875, 0.032125756193742694, 0.26046875, 0.019213329273115587, 0.26046875, 0.006394148914942832 ], "frequencies": [0, 36000, 72000, 108000, 144000, 180000, 216000, 252000, 288000, 324000, 360000, 396000, 432000, 468000, 504000, 540000, 576000, 612000, 648000, 684000, 720000, 756000, 792000, 828000, 864000, 900000, 936000, 972000, 1008000, 1044000, 1080000, 1116000, 1152000, 1188000, 1224000, 1260000, 1296000, 1332000, 1368000, 1404000, 1440000, 1476000, 1512000, 1548000, 1584000, 1620000, 1656000, 1692000, 1728000, 1764000, 1800000, 1836000, 1872000, 1908000, 1944000, 1980000, 2016000, 2052000, 2088000, 2124000, 2160000, 2196000, 2232000, 2268000 ] } } } ] } ] }, "magnetic": {"coil": {"functionalDescription": [{"name": "Primary"}, {"name": "Secondary"}, {"name": "Tertiary"} ] } } })";
        json masJson = json::parse(masString);

        if (masJson["inputs"]["designRequirements"]["insulation"]["cti"] == "GROUP_IIIa") {
            masJson["inputs"]["designRequirements"]["insulation"]["cti"] = "GROUP_IIIA";
        }
        if (masJson["inputs"]["designRequirements"]["insulation"]["cti"] == "GROUP_IIIb") {
            masJson["inputs"]["designRequirements"]["insulation"]["cti"] = "GROUP_IIIB";
        }

        OpenMagnetics::InputsWrapper inputs(masJson["inputs"]);
        OpenMagnetics::MasWrapper masMagnetic;
        // inputs.process_waveforms();

        OpenMagnetics::MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 1);

        CHECK(masMagnetics.size() > 0);

        for (auto [masMagnetic, scoring] : masMagnetics) {
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            if (plot) {
                auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
                auto outFile = outputFilePath;
                std::string filename = "MagneticAdviser_MagneticAdviserJsonLV_" + std::to_string(scoring) + ".svg";
                outFile.append(filename);
                OpenMagnetics::Painter painter(outFile);

                painter.paint_core(masMagnetic.get_mutable_magnetic());
                painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
                painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
                painter.export_svg();
            }
        }
    }

    TEST(Test_MagneticAdviser_Random) {
        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->reset();
        srand (time(NULL));

        int count = 10;
        while (count > 0) {
            std::vector<double> turnsRatios;

            std::vector<int64_t> numberTurns;
            std::vector<OpenMagnetics::IsolationSide> isolationSides;
            // for (size_t windingIndex = 0; windingIndex < std::rand() % 2 + 1UL; ++windingIndex)
            for (size_t windingIndex = 0; windingIndex < std::rand() % 4 + 1UL; ++windingIndex)
            {
                numberTurns.push_back(std::rand() % 100 + 1L);
                isolationSides.push_back(OpenMagnetics::get_isolation_side_from_index(static_cast<size_t>(std::rand() % 10 + 1L)));
            }
            for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
                turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
            }

            double frequency = std::rand() % 1000000 + 10000L;
            double magnetizingInductance = double(std::rand() % 10000) *  1e-6;
            double temperature = 25;
            double peakToPeak = std::rand() % 30 + 1L;
            double dutyCycle = double(std::rand() % 99 + 1L) / 100;
            double dcCurrent = 0;
            if (numberTurns.size() == 1) {
                dcCurrent = std::rand() % 30;
            }

            int waveShapeIndex = std::rand() % magic_enum::enum_count<OpenMagnetics::WaveformLabel>();

            OpenMagnetics::WaveformLabel waveShape = magic_enum::enum_cast<OpenMagnetics::WaveformLabel>(waveShapeIndex).value();
            if (waveShape == OpenMagnetics::WaveformLabel::BIPOLAR_RECTANGULAR ||
                waveShape == OpenMagnetics::WaveformLabel::CUSTOM ||
                waveShape == OpenMagnetics::WaveformLabel::RECTANGULAR ||
                waveShape == OpenMagnetics::WaveformLabel::UNIPOLAR_RECTANGULAR) {
                waveShape = OpenMagnetics::WaveformLabel::TRIANGULAR;
            }

            auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  waveShape,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  dcCurrent,
                                                                                                  turnsRatios);

            inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
            inputs.process_waveforms();

            try {
                OpenMagnetics::MagneticAdviser MagneticAdviser;
                auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 1);
                count--;
            }
            catch (...) {
                for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
                    std::cout << "numberTurns: " << numberTurns[windingIndex] << std::endl;
                }
                for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
                    std::cout << "isolationSides: " << magic_enum::enum_name(isolationSides[windingIndex]) << std::endl;
                }
                std::cout << "frequency: " << frequency << std::endl;
                std::cout << "peakToPeak: " << peakToPeak << std::endl;
                std::cout << "magnetizingInductance: " << magnetizingInductance << std::endl;
                std::cout << "dutyCycle: " << dutyCycle << std::endl;
                std::cout << "dcCurrent: " << dcCurrent << std::endl;
                std::cout << "waveShape: " << magic_enum::enum_name(waveShape) << std::endl;
                CHECK(false);
                return;
            }

        }
    }

    TEST(Test_MagneticAdviser_Random_0) {
        srand (time(NULL));

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {40, 75};
        std::vector<OpenMagnetics::IsolationSide> isolationSides = {OpenMagnetics::IsolationSide::SENARY,
                                                                    OpenMagnetics::IsolationSide::SECONDARY};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 25280;
        double peakToPeak = 28;
        double magnetizingInductance = 0.004539;
        double dutyCycle = 0.68;
        double dcCurrent = 0;
        double temperature = 25;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;

        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process_waveforms();

        OpenMagnetics::MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 1);

        if (masMagnetics.size() > 0) {
            auto masMagnetic = masMagnetics[0].first;
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "MagneticAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagnetic.get_mutable_magnetic());
            painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            painter.export_svg();
        }
    }

    TEST(Test_MagneticAdviser_Random_1) {
        srand (time(NULL));

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {24};
        std::vector<OpenMagnetics::IsolationSide> isolationSides = {OpenMagnetics::IsolationSide::SECONDARY};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 47405;
        double peakToPeak = 25;
        double magnetizingInductance = 0.000831;
        double dutyCycle = 0.05;
        double dcCurrent = 6;
        double temperature = 25;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY;

        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process_waveforms();

        OpenMagnetics::MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 1);

        if (masMagnetics.size() > 0) {
            auto masMagnetic = masMagnetics[0].first;
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "MagneticAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagnetic.get_mutable_magnetic());
            painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            painter.export_svg();
        }
    }

    TEST(Test_MagneticAdviser_Random_2) {
        srand (time(NULL));

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {45, 94};
        std::vector<OpenMagnetics::IsolationSide> isolationSides = {OpenMagnetics::IsolationSide::SECONDARY, OpenMagnetics::IsolationSide::DENARY};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 569910;
        double peakToPeak = 1;
        double magnetizingInductance = 0.00354;
        double dutyCycle = 0.01;
        double dcCurrent = 0;
        double temperature = 25;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::TRIANGULAR;

        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process_waveforms();

        OpenMagnetics::MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 1);

        if (masMagnetics.size() > 0) {
            auto masMagnetic = masMagnetics[0].first;
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "MagneticAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagnetic.get_mutable_magnetic());
            painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            painter.export_svg();
        }
    }

    TEST(Test_MagneticAdviser_Random_3) {
        OpenMagnetics::clear_databases();
        settings->reset();
        srand (time(NULL));

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {78, 47, 1, 100};
        std::vector<OpenMagnetics::IsolationSide> isolationSides = {OpenMagnetics::IsolationSide::QUATERNARY,
                                                                    OpenMagnetics::IsolationSide::TERTIARY,
                                                                    OpenMagnetics::IsolationSide::QUINARY,
                                                                    OpenMagnetics::IsolationSide::UNDENARY};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 811888;
        double peakToPeak = 1;
        double magnetizingInductance = 0.001213;
        double dutyCycle = 0.92;
        double dcCurrent = 0;
        double temperature = 25;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::FLYBACK_SECONDARY;

        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process_waveforms();

        OpenMagnetics::MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 1);

        if (masMagnetics.size() > 0) {
            auto masMagnetic = masMagnetics[0].first;
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "MagneticAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagnetic.get_mutable_magnetic());
            painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            painter.export_svg();
        }
    }

    TEST(Test_MagneticAdviser_Web_0) {
        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->set_use_only_cores_in_stock(true);

        OpenMagnetics::InputsWrapper inputs = json::parse(R"({"designRequirements": {"insulation": {"altitude": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 2000, "minimum": null, "nominal": null }, "cti": "Group II", "insulationType": "Reinforced", "mainSupplyVoltage": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 400, "minimum": null, "nominal": null }, "overvoltageCategory": "OVC-II", "pollutionDegree": "P2", "standards": ["IEC 62368-1" ] }, "isolationSides": ["primary", "secondary", "tertiary" ], "leakageInductance": null, "magnetizingInductance": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.00086 }, "market": null, "maximumDimensions": null, "maximumWeight": null, "name": "My Design Requirements", "operatingTemperature": null, "strayCapacitance": null, "terminalType": null, "topology": null, "turnsRatios": [{"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 7.5 }, {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 8.2 } ], "wiringTechnology": null }, "operatingPoints": [{"conditions": {"ambientRelativeHumidity": null, "ambientTemperature": 42, "cooling": null, "name": null }, "excitationsPerWinding": [{"current": {"harmonics": {"amplitudes": [0.09185791015625004, 0.16258353477989718, 0.11093320057119167, 0.05899976171244, 0.04019583139358277, 0.03853633055544948, 0.029580315036848454, 0.023864828612195715, 0.023362171457053202, 0.019810038012040544, 0.0170566652173707, 0.01683028086934158, 0.014952573418289579, 0.013328356307439995, 0.013210979514049995, 0.012060826394093609, 0.010985890734441771, 0.010922960598071558, 0.010152587942887216, 0.009386239561695625, 0.009354680084182058, 0.008807055146927478, 0.00823147306674625, 0.008219884664203847, 0.007813963766743062, 0.007364651763531841, 0.0073667958811322095, 0.0070566350544148775, 0.006695304572719339, 0.0067075193147442215, 0.006465197241769, 0.006167673146859083, 0.00618770602447329, 0.005995325574340445, 0.005745586287679112, 0.0057720142724450385, 0.005617615641661268, 0.005404620359418445, 0.005436540226690036, 0.005311865059159139, 0.0051277572003770346, 0.0051646110083933714, 0.005063810137351169, 0.004902847424350087, 0.004944320376138756, 0.004863172861122102, 0.004721062807687317, 0.004767022332322001, 0.0047024461235914425, 0.004575918437866753, 0.004626378193656656, 0.004576114721411069, 0.0044626377542477105, 0.004517737729856059, 0.004480144523660421, 0.004377732372497507, 0.004437730265067445, 0.0044116433970340626, 0.0043187218455655405, 0.004383993198211502, 0.004368636843965324, 0.004283948583203358, 0.004354994624223753, 0.004349924130840414 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] }, "processed": {"acEffectiveFrequency": 730484.0527066954, "average": null, "dutyCycle": 0.33, "effectiveFrequency": 676548.4700904511, "label": "Flyback Primary", "offset": 0, "peak": 0.5468750000000006, "peakToPeak": 0.55, "rms": 0.18408843506305278, "thd": 0.9619249449464272 }, "waveform": {"ancillaryLabel": null, "data": [0, 0.013020833333333332, 0.026041666666666664, 0.0390625, 0.05208333333333333, 0.06510416666666666, 0.078125, 0.09114583333333334, 0.10416666666666666, 0.1171875, 0.13020833333333331, 0.14322916666666666, 0.15625, 0.16927083333333334, 0.18229166666666669, 0.1953125, 0.20833333333333331, 0.22135416666666666, 0.234375, 0.24739583333333334, 0.26041666666666663, 0.2734375, 0.2864583333333333, 0.29947916666666663, 0.3125, 0.3255208333333333, 0.3385416666666667, 0.35156250000000006, 0.3645833333333334, 0.37760416666666674, 0.39062500000000017, 0.4036458333333335, 0.4166666666666669, 0.4296875000000002, 0.4427083333333336, 0.45572916666666696, 0.46875000000000033, 0.4817708333333337, 0.4947916666666671, 0.5078125000000004, 0.5208333333333338, 0.5338541666666672, 0.5468750000000006, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], "numberPeriods": null, "time": [0, 7.8125e-8, 1.5625e-7, 2.3437500000000003e-7, 3.125e-7, 3.90625e-7, 4.6875e-7, 5.468750000000001e-7, 6.25e-7, 7.03125e-7, 7.8125e-7, 8.59375e-7, 9.375e-7, 0.0000010156250000000001, 0.0000010937500000000001, 0.0000011718750000000001, 0.00000125, 0.000001328125, 0.00000140625, 0.000001484375, 0.0000015625, 0.000001640625, 0.00000171875, 0.000001796875, 0.000001875, 0.000001953125, 0.0000020312500000000002, 0.0000021093750000000005, 0.0000021875000000000007, 0.000002265625000000001, 0.000002343750000000001, 0.0000024218750000000013, 0.0000025000000000000015, 0.0000025781250000000017, 0.000002656250000000002, 0.000002734375000000002, 0.0000028125000000000023, 0.0000028906250000000025, 0.0000029687500000000027, 0.000003046875000000003, 0.000003125000000000003, 0.0000032031250000000033, 0.0000032812500000000035, 0.0000033593750000000037, 0.000003437500000000004, 0.000003515625000000004, 0.0000035937500000000043, 0.0000036718750000000045, 0.0000037500000000000048, 0.000003828125000000005, 0.000003906250000000005, 0.000003984375000000005, 0.0000040625000000000056, 0.000004140625000000006, 0.000004218750000000006, 0.000004296875000000006, 0.000004375000000000006, 0.000004453125000000007, 0.000004531250000000007, 0.000004609375000000007, 0.000004687500000000007, 0.000004765625000000007, 0.000004843750000000008, 0.000004921875000000008, 0.000005000000000000008, 0.000005078125000000008, 0.0000051562500000000084, 0.000005234375000000009, 0.000005312500000000009, 0.000005390625000000009, 0.000005468750000000009, 0.0000055468750000000095, 0.00000562500000000001, 0.00000570312500000001, 0.00000578125000000001, 0.00000585937500000001, 0.0000059375000000000105, 0.000006015625000000011, 0.000006093750000000011, 0.000006171875000000011, 0.000006250000000000011, 0.0000063281250000000115, 0.000006406250000000012, 0.000006484375000000012, 0.000006562500000000012, 0.000006640625000000012, 0.0000067187500000000125, 0.000006796875000000013, 0.000006875000000000013, 0.000006953125000000013, 0.000007031250000000013, 0.0000071093750000000136, 0.000007187500000000014, 0.000007265625000000014, 0.000007343750000000014, 0.000007421875000000014, 0.000007500000000000015, 0.000007578125000000015, 0.000007656250000000015, 0.000007734375000000015, 0.000007812500000000015, 0.000007890625000000016, 0.000007968750000000016, 0.000008046875000000016, 0.000008125000000000016, 0.000008203125000000016, 0.000008281250000000017, 0.000008359375000000017, 0.000008437500000000017, 0.000008515625000000017, 0.000008593750000000017, 0.000008671875000000018, 0.000008750000000000018, 0.000008828125000000018, 0.000008906250000000018, 0.000008984375000000018, 0.000009062500000000019, 0.000009140625000000019, 0.000009218750000000019, 0.00000929687500000002, 0.00000937500000000002, 0.00000945312500000002, 0.00000953125000000002, 0.00000960937500000002, 0.00000968750000000002, 0.00000976562500000002, 0.00000984375000000002, 0.000009921875000000021 ] } }, "frequency": 100000, "magneticFieldStrength": null, "magneticFluxDensity": {"harmonics": null, "processed": {"acEffectiveFrequency": null, "average": 0.04886606578419124, "dutyCycle": 0.4999999999999997, "effectiveFrequency": null, "label": "Triangular", "offset": 0.04886606578419125, "peak": 0.0977321315683825, "peakToPeak": 0.0977321315683825, "rms": null, "thd": null }, "waveform": {"ancillaryLabel": null, "data": [0, 0.0977321315683825, 0 ], "numberPeriods": null, "time": [0, 0.000005000000000000008, 0.000010000000000000021 ] } }, "magnetizingCurrent": {"harmonics": {"amplitudes": [0.27499999999999936, 0.22295136852075673, 9.138652676247671e-17, 0.024812208134911966, 4.965544356355217e-17, 0.008961158410905185, 2.6247336669766727e-17, 0.004594138958530391, 2.3314386257480727e-17, 0.0027971326973489034, 1.8417688062564944e-17, 0.0018876225747179983, 2.020259692848873e-17, 0.0013646593770613724, 8.2838674091778e-18, 0.001036695055587334, 1.3877787807814463e-17, 0.0008176648182778223, 1.693968269140533e-17, 0.0006642449307756816, 6.669234702534747e-18, 0.000552698506300629, 1.3315671033185686e-17, 0.00046913630363004693, 6.414434300065929e-18, 0.00040499257758880064, 1.0921946309111329e-17, 0.00035475251819120247, 6.079119163783553e-18, 0.00031473605872487254, 6.936742671150331e-18, 0.00028241198689922184, 0, 0.00025599367330819076, 6.816639972961569e-18, 0.0002341916374268395, 1.745648778780573e-18, 0.0002160570955530902, 3.614121502139124e-18, 0.0002008801898065849, 8.466014217699864e-18, 0.00018812211783339133, 7.694657952990434e-18, 0.00017736886649857841, 8.120827158253018e-18, 0.00016829905656025996, 9.577499028380518e-18, 0.00016066121193102052, 1.3877787807814475e-17, 0.00015425745363482706, 6.591048062130357e-18, 0.00014893165806567773, 8.910189550509116e-19, 0.0001445607749454254, 3.71282312497132e-18, 0.00014104842315542927, 1.0676750772951972e-17, 0.00013832016090481892, 6.646156529573393e-18, 0.00013632001373806244, 7.513753614217868e-18, 0.0001350079724723105, 2.8701438907104713e-18, 0.00013435826395631623 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] }, "processed": {"acEffectiveFrequency": null, "average": 0.2749999999999994, "dutyCycle": 0.4999999999999997, "effectiveFrequency": null, "label": "Triangular", "offset": 0.2749999999999996, "peak": 0.5499999999999992, "peakToPeak": 0.5499999999999992, "rms": null, "thd": null }, "waveform": {"ancillaryLabel": null, "data": [0, 0.00859375, 0.0171875, 0.025781250000000002, 0.034375, 0.04296875, 0.051562500000000004, 0.06015625000000001, 0.06875, 0.07734375, 0.0859375, 0.09453125000000001, 0.10312500000000001, 0.11171875, 0.12031250000000002, 0.12890625, 0.1375, 0.14609375000000002, 0.1546875, 0.16328125000000002, 0.171875, 0.18046875, 0.18906250000000002, 0.19765625, 0.20625000000000002, 0.21484375000000003, 0.2234375, 0.23203125000000005, 0.2406250000000001, 0.24921875000000007, 0.2578125000000001, 0.2664062500000001, 0.27500000000000013, 0.2835937500000002, 0.2921875000000002, 0.3007812500000002, 0.3093750000000003, 0.31796875000000024, 0.32656250000000026, 0.3351562500000003, 0.34375000000000033, 0.3523437500000004, 0.3609375000000004, 0.3695312500000004, 0.37812500000000043, 0.38671875000000044, 0.39531250000000046, 0.4039062500000005, 0.41250000000000053, 0.42109375000000054, 0.4296875000000006, 0.43828125000000057, 0.4468750000000006, 0.45546875000000064, 0.46406250000000066, 0.4726562500000007, 0.48125000000000073, 0.4898437500000007, 0.49843750000000075, 0.5070312500000008, 0.5156250000000008, 0.5242187500000008, 0.5328125000000008, 0.5414062500000009, 0.5499999999999992, 0.5414062499999992, 0.5328124999999991, 0.5242187499999992, 0.5156249999999991, 0.507031249999999, 0.4984374999999991, 0.4898437499999991, 0.48124999999999907, 0.472656249999999, 0.464062499999999, 0.4554687499999989, 0.4468749999999989, 0.43828124999999896, 0.4296874999999989, 0.4210937499999989, 0.4124999999999988, 0.4039062499999988, 0.3953124999999988, 0.38671874999999883, 0.37812499999999877, 0.36953124999999876, 0.3609374999999987, 0.3523437499999987, 0.34374999999999867, 0.33515624999999866, 0.32656249999999865, 0.3179687499999986, 0.30937499999999857, 0.30078124999999856, 0.29218749999999855, 0.2835937499999985, 0.27499999999999847, 0.26640624999999846, 0.25781249999999845, 0.24921874999999843, 0.24062499999999837, 0.23203124999999836, 0.22343749999999837, 0.2148437499999983, 0.2062499999999983, 0.1976562499999983, 0.18906249999999825, 0.18046874999999823, 0.17187499999999825, 0.16328124999999818, 0.15468749999999817, 0.1460937499999981, 0.13749999999999812, 0.1289062499999981, 0.12031249999999806, 0.11171874999999805, 0.10312499999999805, 0.094531249999998, 0.08593749999999799, 0.07734374999999799, 0.06874999999999792, 0.06015624999999793, 0.05156249999999793, 0.04296874999999786, 0.034374999999997866, 0.025781249999997802, 0.017187499999997802, 0.008593749999997803 ], "numberPeriods": null, "time": [0, 7.8125e-8, 1.5625e-7, 2.3437500000000003e-7, 3.125e-7, 3.90625e-7, 4.6875e-7, 5.468750000000001e-7, 6.25e-7, 7.03125e-7, 7.8125e-7, 8.59375e-7, 9.375e-7, 0.0000010156250000000001, 0.0000010937500000000001, 0.0000011718750000000001, 0.00000125, 0.000001328125, 0.00000140625, 0.000001484375, 0.0000015625, 0.000001640625, 0.00000171875, 0.000001796875, 0.000001875, 0.000001953125, 0.0000020312500000000002, 0.0000021093750000000005, 0.0000021875000000000007, 0.000002265625000000001, 0.000002343750000000001, 0.0000024218750000000013, 0.0000025000000000000015, 0.0000025781250000000017, 0.000002656250000000002, 0.000002734375000000002, 0.0000028125000000000023, 0.0000028906250000000025, 0.0000029687500000000027, 0.000003046875000000003, 0.000003125000000000003, 0.0000032031250000000033, 0.0000032812500000000035, 0.0000033593750000000037, 0.000003437500000000004, 0.000003515625000000004, 0.0000035937500000000043, 0.0000036718750000000045, 0.0000037500000000000048, 0.000003828125000000005, 0.000003906250000000005, 0.000003984375000000005, 0.0000040625000000000056, 0.000004140625000000006, 0.000004218750000000006, 0.000004296875000000006, 0.000004375000000000006, 0.000004453125000000007, 0.000004531250000000007, 0.000004609375000000007, 0.000004687500000000007, 0.000004765625000000007, 0.000004843750000000008, 0.000004921875000000008, 0.000005000000000000008, 0.000005078125000000008, 0.0000051562500000000084, 0.000005234375000000009, 0.000005312500000000009, 0.000005390625000000009, 0.000005468750000000009, 0.0000055468750000000095, 0.00000562500000000001, 0.00000570312500000001, 0.00000578125000000001, 0.00000585937500000001, 0.0000059375000000000105, 0.000006015625000000011, 0.000006093750000000011, 0.000006171875000000011, 0.000006250000000000011, 0.0000063281250000000115, 0.000006406250000000012, 0.000006484375000000012, 0.000006562500000000012, 0.000006640625000000012, 0.0000067187500000000125, 0.000006796875000000013, 0.000006875000000000013, 0.000006953125000000013, 0.000007031250000000013, 0.0000071093750000000136, 0.000007187500000000014, 0.000007265625000000014, 0.000007343750000000014, 0.000007421875000000014, 0.000007500000000000015, 0.000007578125000000015, 0.000007656250000000015, 0.000007734375000000015, 0.000007812500000000015, 0.000007890625000000016, 0.000007968750000000016, 0.000008046875000000016, 0.000008125000000000016, 0.000008203125000000016, 0.000008281250000000017, 0.000008359375000000017, 0.000008437500000000017, 0.000008515625000000017, 0.000008593750000000017, 0.000008671875000000018, 0.000008750000000000018, 0.000008828125000000018, 0.000008906250000000018, 0.000008984375000000018, 0.000009062500000000019, 0.000009140625000000019, 0.000009218750000000019, 0.00000929687500000002, 0.00000937500000000002, 0.00000945312500000002, 0.00000953125000000002, 0.00000960937500000002, 0.00000968750000000002, 0.00000976562500000002, 0.00000984375000000002, 0.000009921875000000021 ] } }, "name": "Primary winding excitation", "voltage": {"harmonics": {"amplitudes": [0.2687500000000034, 78.27459402394327, 40.25330713184367, 1.4938054445226292, 18.99816174625034, 16.53907513455175, 1.4960597290909847, 10.52193687000004, 10.60588337564493, 1.4998323370832707, 7.124948744925913, 7.906132484799749, 1.5051467332495563, 5.290144683108919, 6.36156775076611, 1.512036305031938, 4.138210406040141, 5.36091892012817, 1.5205449194517633, 3.345138664057968, 4.659663922070232, 1.5307276689596565, 2.763604439072038, 4.1408468806075795, 1.5426518307315582, 2.316949035987901, 3.741483514643341, 1.556398071861931, 1.9613162141323381, 3.424645836484155, 1.572061942683373, 1.6697726425283888, 3.1672491240647447, 1.5897557126829684, 1.4248338912738232, 2.954119814429578, 1.6096106190467687, 1.2146249767516817, 2.7748732925030497, 1.6317796179173467, 1.0307669990419286, 2.6221623682054127, 1.6564407546007427, 0.8671465135513805, 2.4906422561346835, 1.6838013034499004, 0.7191641870611972, 2.3763319237335634, 1.7141028741352993, 0.5832567645981096, 2.2762053377540736, 1.7476277429945672, 0.45658126987367614, 2.1879215212119307, 1.7847067526353455, 0.3367986140087129, 2.109641370402089, 1.8257292394763993, 0.22191949733106117, 2.0399003544286045, 1.8711556115945902, 0.11018970942035004, 1.977518172190299, 1.9215334294233297 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] }, "processed": {"acEffectiveFrequency": 638492.1664714317, "average": -0.26875000000005356, "dutyCycle": 0.32812499999999967, "effectiveFrequency": 638489.6210029337, "label": "Rectangular", "offset": 0, "peak": 96.03333333333332, "peakToPeak": 143.33333333333331, "rms": 67.29992261015853, "thd": 0.6917104252878409 }, "waveform": {"ancillaryLabel": null, "data": [-47.3, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, 96.03333333333332, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3, -47.3 ], "numberPeriods": null, "time": [0, 7.8125e-8, 1.5625e-7, 2.3437500000000003e-7, 3.125e-7, 3.90625e-7, 4.6875e-7, 5.468750000000001e-7, 6.25e-7, 7.03125e-7, 7.8125e-7, 8.59375e-7, 9.375e-7, 0.0000010156250000000001, 0.0000010937500000000001, 0.0000011718750000000001, 0.00000125, 0.000001328125, 0.00000140625, 0.000001484375, 0.0000015625, 0.000001640625, 0.00000171875, 0.000001796875, 0.000001875, 0.000001953125, 0.0000020312500000000002, 0.0000021093750000000005, 0.0000021875000000000007, 0.000002265625000000001, 0.000002343750000000001, 0.0000024218750000000013, 0.0000025000000000000015, 0.0000025781250000000017, 0.000002656250000000002, 0.000002734375000000002, 0.0000028125000000000023, 0.0000028906250000000025, 0.0000029687500000000027, 0.000003046875000000003, 0.000003125000000000003, 0.0000032031250000000033, 0.0000032812500000000035, 0.0000033593750000000037, 0.000003437500000000004, 0.000003515625000000004, 0.0000035937500000000043, 0.0000036718750000000045, 0.0000037500000000000048, 0.000003828125000000005, 0.000003906250000000005, 0.000003984375000000005, 0.0000040625000000000056, 0.000004140625000000006, 0.000004218750000000006, 0.000004296875000000006, 0.000004375000000000006, 0.000004453125000000007, 0.000004531250000000007, 0.000004609375000000007, 0.000004687500000000007, 0.000004765625000000007, 0.000004843750000000008, 0.000004921875000000008, 0.000005000000000000008, 0.000005078125000000008, 0.0000051562500000000084, 0.000005234375000000009, 0.000005312500000000009, 0.000005390625000000009, 0.000005468750000000009, 0.0000055468750000000095, 0.00000562500000000001, 0.00000570312500000001, 0.00000578125000000001, 0.00000585937500000001, 0.0000059375000000000105, 0.000006015625000000011, 0.000006093750000000011, 0.000006171875000000011, 0.000006250000000000011, 0.0000063281250000000115, 0.000006406250000000012, 0.000006484375000000012, 0.000006562500000000012, 0.000006640625000000012, 0.0000067187500000000125, 0.000006796875000000013, 0.000006875000000000013, 0.000006953125000000013, 0.000007031250000000013, 0.0000071093750000000136, 0.000007187500000000014, 0.000007265625000000014, 0.000007343750000000014, 0.000007421875000000014, 0.000007500000000000015, 0.000007578125000000015, 0.000007656250000000015, 0.000007734375000000015, 0.000007812500000000015, 0.000007890625000000016, 0.000007968750000000016, 0.000008046875000000016, 0.000008125000000000016, 0.000008203125000000016, 0.000008281250000000017, 0.000008359375000000017, 0.000008437500000000017, 0.000008515625000000017, 0.000008593750000000017, 0.000008671875000000018, 0.000008750000000000018, 0.000008828125000000018, 0.000008906250000000018, 0.000008984375000000018, 0.000009062500000000019, 0.000009140625000000019, 0.000009218750000000019, 0.00000929687500000002, 0.00000937500000000002, 0.00000945312500000002, 0.00000953125000000002, 0.00000960937500000002, 0.00000968750000000002, 0.00000976562500000002, 0.00000984375000000002, 0.000009921875000000021 ] } } }, {"current": {"harmonics": {"amplitudes": [1.331841767723876, 1.5978437326305774, 0.5814203703679136, 0.422682846269918, 0.3357339897017588, 0.24494681037000157, 0.21191462908744038, 0.18801596328155598, 0.15588992203480992, 0.1419165456125347, 0.1310820681222257, 0.11478576298308843, 0.10711467401070128, 0.10104516981227131, 0.09122305193126073, 0.08639547837402142, 0.08257801501466873, 0.0760250641959211, 0.07272221550131976, 0.07014374133895233, 0.06546938790559159, 0.06308015385360839, 0.061255812980550245, 0.05776025934348251, 0.05596289609270588, 0.05463232568083056, 0.05192540562509951, 0.05053450484458919, 0.04954587549843613, 0.04739307515142994, 0.04629437227456433, 0.045553528464656295, 0.04380555890610046, 0.04292490892933754, 0.04237065959121802, 0.04092809879015255, 0.04021541120437753, 0.03980654237733284, 0.03860074516438956, 0.03802109284385031, 0.03772910362032111, 0.03671127625211673, 0.0362397040479845, 0.03604453630130419, 0.03517917878789851, 0.03479753740343691, 0.034684963952665966, 0.03394577636445149, 0.03364072217763382, 0.0336007002714556, 0.03296794337692874, 0.03272963965589644, 0.03275525215925737, 0.032214001245385364, 0.03203525860776655, 0.03212202911541702, 0.03166099669766672, 0.03153669973459098, 0.03168215731738724, 0.03129289134787348, 0.03121961986783353, 0.031423040555606546, 0.03109938021000907, 0.0310751699340615 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] }, "processed": {"acEffectiveFrequency": 631602.8673449805, "average": null, "dutyCycle": 0.33, "effectiveFrequency": 515359.4457714902, "label": "Flyback Secondary", "offset": 0, "peak": 3.964552238805968, "peakToPeak": 4, "rms": 1.8817041929560256, "thd": 0.6196863660749166 }, "waveform": {"ancillaryLabel": null, "data": [0, 0, 4, 0, 0 ], "numberPeriods": null, "time": [0, 0.0000033000000000000006, 0.0000033000000000000006, 0.00001, 0.00001 ] } }, "frequency": 100000, "magneticFieldStrength": null, "magneticFluxDensity": null, "magnetizingCurrent": {"harmonics": {"amplitudes": [1.9999999999999953, 1.621464498332776, 7.038312373140319e-16, 0.1804524227993598, 3.662638187907386e-16, 0.06517206117021954, 1.9182595835467152e-16, 0.03341191969840281, 1.967149070383611e-16, 0.02034278325344657, 1.350828826864373e-16, 0.013728164179767275, 2.1072841793042454e-16, 0.009924795469537248, 1.0907257763864025e-16, 0.00753960040427154, 2.476288615702221e-16, 0.005946653223838669, 9.099075076946043e-17, 0.004830872223823178, 6.790339492243881e-17, 0.004019625500368258, 9.693022546373344e-17, 0.003411900390036779, 3.802880834339934e-18, 0.0029454005642822224, 6.623718649805912e-17, 0.002580018314117822, 3.9074282406338975e-17, 0.002288989517999055, 3.6696069249880155e-17, 0.0020539053592668333, 1.570092458683776e-16, 0.001861772169514179, 2.4159946427735383e-17, 0.0017032119085588196, 6.644505947526663e-17, 0.0015713243312952172, 5.968027579395866e-17, 0.001460946834956927, 6.422169018553131e-17, 0.0013681608569701188, 4.9392738563415026e-17, 0.0012899553927169048, 3.740986020933547e-17, 0.0012239931386200823, 3.490090048248212e-17, 0.0011684451776802197, 1.7598441346534375e-17, 0.001121872390071435, 4.0397251965478305e-17, 0.0010831393313867413, 3.2007004974977725e-17, 0.0010513510905121957, 1.871613877384591e-17, 0.0010258067138577028, 1.7469651119485624e-17, 0.001005964806580533, 3.047165160071233e-17, 0.0009914182817314099, 8.499268590616167e-17, 0.0009818761634349943, 3.696179851046128e-17, 0.0009771510105917036 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] }, "processed": {"acEffectiveFrequency": 110746.40291779555, "average": 1.9999999999999976, "dutyCycle": 0.4999999999999997, "effectiveFrequency": 70052.43250784607, "label": "Triangular", "offset": 1.999999999999997, "peak": 3.999999999999994, "peakToPeak": 3.999999999999994, "rms": 2.3095420271127347, "thd": 0.12151487440704967 }, "waveform": {"ancillaryLabel": null, "data": [0, 3.999999999999994, 0 ], "numberPeriods": null, "time": [0, 0.000005000000000000008, 0.000010000000000000021 ] } }, "name": "Primary winding excitation", "voltage": {"harmonics": {"amplitudes": [0.022500000000000256, 6.553221825260367, 3.3700443180148194, 0.12506278140189475, 1.5905437741046795, 1.3846667554508447, 0.12525151220296624, 0.8809063426046547, 0.8879344221470176, 0.1255673584534832, 0.5965073367844952, 0.6619087661692814, 0.12601228464414885, 0.4428958339347003, 0.5325963698315814, 0.12658908600267368, 0.3464548246917328, 0.4488211188944513, 0.12730143511689168, 0.2800581207118298, 0.3901113981268103, 0.1281539443780178, 0.23137153443393804, 0.3466755527950532, 0.12915224629380498, 0.19397712859433597, 0.3132404802957216, 0.13030309438844082, 0.1642032179273585, 0.2867145351475106, 0.1316144882246547, 0.13979491890935458, 0.2651650429449554, 0.1330958271083416, 0.1192884188043201, 0.24732165888247626, 0.13475809833879931, 0.10168953293734996, 0.2323149733258367, 0.13661410754656847, 0.08629677201281269, 0.219529872686965, 0.13867876085029474, 0.07259831276244126, 0.20851888656011303, 0.1409694114516196, 0.06020909473070487, 0.1989487191962986, 0.1435062871369088, 0.048830798896586015, 0.19056602827708508, 0.14631302034373123, 0.03822540864058689, 0.18317482503169658, 0.14941730952295895, 0.028197093265845675, 0.17662113798715165, 0.15285175028174516, 0.018579306753298146, 0.17078235525448782, 0.1566548884125703, 0.009225184974726776, 0.16555966092756036, 0.1608725661842781 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] }, "processed": {"acEffectiveFrequency": 638492.1664714315, "average": -0.022499999999998, "dutyCycle": 0.32812499999999967, "effectiveFrequency": 638489.6210029334, "label": "Rectangular", "offset": 0, "peak": 8.04, "peakToPeak": 12, "rms": 5.63441212550164, "thd": 0.6917104252878407 }, "waveform": {"ancillaryLabel": null, "data": [-3.96, 8.04, 8.04, -3.96, -3.96 ], "numberPeriods": null, "time": [0, 0, 0.0000033000000000000006, 0.0000033000000000000006, 0.00001 ] } } }, {"current": {"harmonics": {"amplitudes": [0.332960441930969, 0.39946093315764436, 0.1453550925919784, 0.1056707115674795, 0.0839334974254397, 0.06123670259250039, 0.052978657271860095, 0.047003990820388995, 0.03897248050870248, 0.035479136403133674, 0.032770517030556424, 0.028696440745772107, 0.02677866850267532, 0.025261292453067827, 0.022805762982815184, 0.021598869593505356, 0.020644503753667182, 0.019006266048980276, 0.01818055387532994, 0.017535935334738083, 0.016367346976397898, 0.015770038463402098, 0.015313953245137561, 0.014440064835870627, 0.01399072402317647, 0.01365808142020764, 0.012981351406274877, 0.012633626211147298, 0.012386468874609032, 0.011848268787857484, 0.011573593068641082, 0.011388382116164074, 0.010951389726525115, 0.010731227232334386, 0.010592664897804504, 0.010232024697538138, 0.010053852801094382, 0.00995163559433321, 0.00965018629109739, 0.009505273210962577, 0.009432275905080277, 0.009177819063029183, 0.009059926011996126, 0.009011134075326048, 0.008794794696974627, 0.008699384350859228, 0.008671240988166492, 0.008486444091112872, 0.008410180544408456, 0.0084001750678639, 0.008241985844232184, 0.00818240991397411, 0.008188813039814342, 0.008053500311346341, 0.008008814651941637, 0.008030507278854255, 0.00791524917441668, 0.007884174933647745, 0.00792053932934681, 0.00782322283696837, 0.007804904966958383, 0.007855760138901637, 0.007774845052502267, 0.007768792483515375 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] }, "processed": {"acEffectiveFrequency": 631602.8673449805, "average": null, "dutyCycle": 0.33, "effectiveFrequency": 515359.4457714902, "label": "Flyback Secondary", "offset": 0, "peak": 0.991138059701492, "peakToPeak": 1, "rms": 0.4704260482390064, "thd": 0.6196863660749166 }, "waveform": {"ancillaryLabel": null, "data": [0, 0, 1, 0, 0 ], "numberPeriods": null, "time": [0, 0.0000033000000000000006, 0.0000033000000000000006, 0.00001, 0.00001 ] } }, "frequency": 100000, "magneticFieldStrength": null, "magneticFluxDensity": null, "magnetizingCurrent": {"harmonics": {"amplitudes": [0.49999999999999883, 0.405366124583194, 1.7595780932850797e-16, 0.04511310569983995, 9.156595469768465e-17, 0.016293015292554884, 4.795648958866788e-17, 0.008352979924600703, 4.917872675959028e-17, 0.005085695813361642, 3.3770720671609327e-17, 0.0034320410449418188, 5.2682104482606134e-17, 0.002481198867384312, 2.7268144409660063e-17, 0.001884900101067885, 6.190721539255553e-17, 0.0014866633059596672, 2.2747687692365106e-17, 0.0012077180559557945, 1.6975848730609703e-17, 0.0010049063750920645, 2.423255636593336e-17, 0.0008529750975091948, 9.507202085849835e-19, 0.0007363501410705556, 1.655929662451478e-17, 0.0006450045785294555, 9.768570601584744e-18, 0.0005722473794997638, 9.174017312470039e-18, 0.0005134763398167083, 3.92523114670944e-17, 0.00046544304237854475, 6.039986606933846e-18, 0.0004258029771397049, 1.6611264868816657e-17, 0.0003928310828238043, 1.4920068948489666e-17, 0.00036523670873923176, 1.6055422546382828e-17, 0.0003420402142425297, 1.2348184640853756e-17, 0.0003224888481792262, 9.352465052333868e-18, 0.00030599828465502056, 8.72522512062053e-18, 0.0002921112944200549, 4.399610336633594e-18, 0.0002804680975178588, 1.0099312991369576e-17, 0.0002707848328466853, 8.001751243744431e-18, 0.0002628377726280489, 4.6790346934614775e-18, 0.0002564516784644257, 4.367412779871406e-18, 0.00025149120164513327, 7.617912900178082e-18, 0.00024785457043285247, 2.1248171476540418e-17, 0.00024546904085874857, 9.24044962761532e-18, 0.0002442877526479259 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] }, "processed": {"acEffectiveFrequency": 110746.40291779555, "average": 0.4999999999999994, "dutyCycle": 0.4999999999999997, "effectiveFrequency": 70052.43250784607, "label": "Triangular", "offset": 0.4999999999999992, "peak": 0.9999999999999984, "peakToPeak": 0.9999999999999984, "rms": 0.5773855067781837, "thd": 0.12151487440704967 }, "waveform": {"ancillaryLabel": null, "data": [0, 0.9999999999999984, 0 ], "numberPeriods": null, "time": [0, 0.000005000000000000008, 0.000010000000000000021 ] } }, "name": "Primary winding excitation", "voltage": {"harmonics": {"amplitudes": [0.022500000000000256, 6.553221825260367, 3.3700443180148194, 0.12506278140189475, 1.5905437741046795, 1.3846667554508447, 0.12525151220296624, 0.8809063426046547, 0.8879344221470176, 0.1255673584534832, 0.5965073367844952, 0.6619087661692814, 0.12601228464414885, 0.4428958339347003, 0.5325963698315814, 0.12658908600267368, 0.3464548246917328, 0.4488211188944513, 0.12730143511689168, 0.2800581207118298, 0.3901113981268103, 0.1281539443780178, 0.23137153443393804, 0.3466755527950532, 0.12915224629380498, 0.19397712859433597, 0.3132404802957216, 0.13030309438844082, 0.1642032179273585, 0.2867145351475106, 0.1316144882246547, 0.13979491890935458, 0.2651650429449554, 0.1330958271083416, 0.1192884188043201, 0.24732165888247626, 0.13475809833879931, 0.10168953293734996, 0.2323149733258367, 0.13661410754656847, 0.08629677201281269, 0.219529872686965, 0.13867876085029474, 0.07259831276244126, 0.20851888656011303, 0.1409694114516196, 0.06020909473070487, 0.1989487191962986, 0.1435062871369088, 0.048830798896586015, 0.19056602827708508, 0.14631302034373123, 0.03822540864058689, 0.18317482503169658, 0.14941730952295895, 0.028197093265845675, 0.17662113798715165, 0.15285175028174516, 0.018579306753298146, 0.17078235525448782, 0.1566548884125703, 0.009225184974726776, 0.16555966092756036, 0.1608725661842781 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] }, "processed": {"acEffectiveFrequency": 638492.1664714315, "average": -0.022499999999998, "dutyCycle": 0.32812499999999967, "effectiveFrequency": 638489.6210029334, "label": "Rectangular", "offset": 0, "peak": 8.04, "peakToPeak": 12, "rms": 5.63441212550164, "thd": 0.6917104252878407 }, "waveform": {"ancillaryLabel": null, "data": [-3.96, 8.04, 8.04, -3.96, -3.96 ], "numberPeriods": null, "time": [0, 0, 0.0000033000000000000006, 0.0000033000000000000006, 0.00001 ] } } } ], "name": "Operating Point No. 1" } ] })");

        OpenMagnetics::MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 2);

        for (auto [mas, scoring] : masMagnetics) {
            std::string name = mas.get_magnetic().get_manufacturer_info().value().get_reference().value();
            auto masMagnetic = mas;
            auto result = OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            name = std::filesystem::path(std::regex_replace(std::string(name), std::regex(" "), "_")).string();
            name = std::filesystem::path(std::regex_replace(std::string(name), std::regex("\\,"), "_")).string();
            name = std::filesystem::path(std::regex_replace(std::string(name), std::regex("\\."), "_")).string();
            name = std::filesystem::path(std::regex_replace(std::string(name), std::regex("\\:"), "_")).string();
            name = std::filesystem::path(std::regex_replace(std::string(name), std::regex("\\/"), "_")).string();
            std::string filename = "MagneticAdviser_" + name + ".svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagnetic.get_mutable_magnetic());
            painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            painter.export_svg();
        }
    }

    TEST(Test_MagneticAdviser_Web_1) {
        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->set_use_only_cores_in_stock(true);

        OpenMagnetics::InputsWrapper inputs = json::parse(R"({"designRequirements": {"insulation": {"altitude": {"maximum": 2000 }, "cti": "Group I", "pollutionDegree": "P2", "overvoltageCategory": "OVC-II", "insulationType": "Reinforced", "mainSupplyVoltage": {"maximum": 400 }, "standards": ["IEC 62368-1" ] }, "magnetizingInductance": {"nominal": 0.00039999999999999996 }, "name": "My Design Requirements", "turnsRatios": [{"nominal": 10 } ] }, "operatingPoints": [{"name": "Operating Point No. 1", "conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"ancillaryLabel": null, "data": [-2.5, 2.5, -2.5 ], "numberPeriods": null, "time": [0, 0.0000025, 0.00001 ] }, "processed": {"dutyCycle": 0.25, "peakToPeak": 5, "offset": 0, "label": "Triangular", "acEffectiveFrequency": 128062.87364573497, "effectiveFrequency": 128062.87364573497, "peak": 2.499999999999999, "rms": 1.4438454453780356, "thd": 0.37655984422983674 }, "harmonics": {"amplitudes": [5.289171878253285e-15, 1.9109142370405827, 0.6760173538930697, 0.21266521973828106, 1.0225829328048648e-15, 0.07680601066227859, 0.07559762456783126, 0.039376324985334776, 4.573645088368781e-16, 0.023974199977867476, 0.027568116438312952, 0.01617879664125947, 2.971250334652087e-16, 0.011696483630665971, 0.014340785210409782, 0.008885504288828753, 2.2404804069938507e-16, 0.0070081980332356445, 0.008903579153043517, 0.0056932375141851765, 2.0060103614117677e-16, 0.004737174081901557, 0.006158131362616469, 0.004020963170880106, 1.6456472022598636e-16, 0.0034711878538580374, 0.004586639500385065, 0.0030405807424970697, 8.639502758732636e-17, 0.0026976000170703556, 0.003608946997642845, 0.00242055067908846, 1.0574701291099621e-16, 0.0021941195434799954, 0.0029646243339663376, 0.00200725448389949, 1.078627841406342e-16, 0.001851823483504021, 0.002522860567684626, 0.0017217423565853373, 8.393610561532662e-17, 0.0016123930328628207, 0.002212325372118857, 0.0015202270093639689, 1.1562789535401302e-16, 0.0014424897474069642, 0.0019916921862528536, 0.0013770258476373952, 7.108157805209173e-17, 0.0013221392910760783, 0.0018359780830164152, 0.0012764919436558824, 9.44838447818231e-17, 0.0012390291425152763, 0.0017297261022555157, 0.0012089248059259444, 1.6792517082756795e-16, 0.0011855408939470595, 0.001663417261065822, 0.0011683976500078048, 8.640837577205282e-17, 0.0011571521557507119, 0.0016315323036584783, 0.0011515835097210575 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-533.3333333333333, 1599.9999999999998, 1599.9999999999998, -533.3333333333333, -533.3333333333333 ], "numberPeriods": null, "time": [0, 0, 0.0000025, 0.0000025, 0.00001 ] }, "processed": {"acEffectiveFrequency": 690362.0785055002, "average": -3.907985046680551e-13, "dutyCycle": 0.25, "effectiveFrequency": 690304.6703245766, "label": "Rectangular", "offset": 533.3333333333334, "peak": 1599.9999999999998, "peakToPeak": 2133.333333333333, "rms": 914.0872800534738, "thd": 0.9507084995254543 }, "harmonics": {"amplitudes": [16.666666666666664, 936.5743366559599, 678.5155874995734, 343.10398442523496, 33.33333333333333, 167.53211803338075, 224.71508018049997, 159.4087857982368, 33.33333333333333, 81.3927569897301, 133.0741261256697, 108.74207884995286, 33.33333333333333, 47.777866517707096, 93.16042574968293, 84.67373936814097, 33.33333333333333, 29.60335188747621, 70.47741191828814, 70.39622283854739, 33.33333333333333, 18.037041704179572, 55.61330685278376, 60.79168042183305, 33.33333333333333, 9.896923507307392, 44.944797116224166, 53.772679523931075, 33.33333333333333, 3.7542162789560147, 36.77766585778291, 48.326808642191196, 33.33333333333333, 1.1295053442673506, 30.211572300638238, 43.9020428459011, 33.33333333333333, 5.175841098510636, 24.721684875734596, 40.17025212715429, 33.33333333333333, 8.644545010278817, 19.979231122730887, 36.92259416373016, 33.33333333333333, 11.70597127653338, 15.765492529710798, 34.018188678401664, 33.33333333333333, 14.478186559464497, 11.926857377150895, 31.356776956938962, 33.33333333333333, 17.047466028773496, 8.349565339710303, 28.863097293035203, 33.33333333333333, 19.480403852519398, 4.944532917945068, 26.477336144614267, 33.33333333333333, 21.831581714888173, 1.6375616589827473, 24.148842611619322 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] } } }, {"current": {"harmonics": {"amplitudes": [0.0000010000000539855947, 19.037231223396237, 6.897984052452256, 2.326563601915805, 0.16964245489488294, 0.6827105539690119, 0.7639584293115499, 0.46512617224180786, 0.08523164109549347, 0.1853099567245175, 0.27319994312359214, 0.2033816061551562, 0.057281247031537384, 0.07597930993201436, 0.13794353928011985, 0.11736286786731481, 0.04345071383027455, 0.036440556746865825, 0.08222723635075783, 0.07810340034727467, 0.03527361775824142, 0.018254330742417037, 0.053961824922667094, 0.056682647846879576, 0.02992937223599544, 0.008573935375330551, 0.03763658083416089, 0.0435994370579927, 0.026210662569045903, 0.002867140660712981, 0.027318562880531572, 0.03496118462997451, 0.023515356873512087, 0.0007779637595615999, 0.020339562018896258, 0.02892190405184444, 0.021510534863982397, 0.0032750781058732345, 0.015353525869835107, 0.024509224663769844, 0.01999816717442588, 0.005104092422150859, 0.011619480338521122, 0.02116836728711663, 0.01885414149705068, 0.006537395429160153, 0.008699656143415562, 0.018561931478014347, 0.017997874963290205, 0.007740940114641076, 0.006318925086584056, 0.016473329652684467, 0.017376076319907764, 0.008823512955740787, 0.004293742261433124, 0.014756505062970141, 0.016953627505956508, 0.00986273033802779, 0.002493499634464789, 0.013307931595332303, 0.016708323451549476, 0.010919968097470808, 0.0008178599364280956, 0.01204994245895335 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] }, "processed": {"acEffectiveFrequency": 129404.9130926142, "average": -5.3734794391857577e-14, "dutyCycle": 0.2421874999999996, "effectiveFrequency": 129404.91309261405, "label": "Triangular", "offset": -0.000001, "peak": 24.999998999999956, "peakToPeak": 49.999999999999986, "rms": 14.438555983981303, "thd": 0.3878858672873039 }, "waveform": {"ancillaryLabel": null, "data": [-25.000000999999994, 24.999998999999992, -25.000000999999994 ], "numberPeriods": null, "time": [0, 0.000002421874999999996, 0.00001 ] } }, "frequency": 100000, "magneticFieldStrength": null, "magneticFluxDensity": null, "magnetizingCurrent": null, "name": "Primary winding excitation", "voltage": {"harmonics": {"amplitudes": [9.16666666666667, 114.01883606805356, 84.50783707683871, 45.4933446767725, 8.293206055601647, 17.4992480352357, 27.17398286677801, 22.031939878495663, 8.173210670026926, 6.406662607033046, 15.123340250417916, 15.15411908517724, 7.974502797768419, 1.9490359117489253, 9.560623024349578, 11.563504110131417, 7.698996104260727, 0.504510858012843, 6.182373666688366, 9.166959370632709, 7.349343869569635, 2.0540091893695203, 3.8205457805728327, 7.332930326725261, 6.928913435854551, 3.093876715391758, 2.0304195076924505, 5.808981955821909, 6.441753778022821, 3.8000361805623823, 0.6081443012940038, 4.4773238922780285, 5.892556509887897, 4.262766710863713, 0.5511898658329515, 3.2779337521433516, 5.286610701363717, 4.53337401311496, 1.5058617370415037, 2.178953825099018, 4.629751941830022, 4.643526917515465, 2.2899470149511534, 1.1635818333216585, 3.928306140216652, 4.614259582955264, 2.924044975741664, 0.22363428278401298, 3.189028603042417, 4.46053575052874, 3.420845617443641, 0.643874494083472, 2.4190389771205214, 4.1937095470147, 3.788199527266002, 1.4386641980705313, 1.6257526834677354, 3.822904217686624, 4.03087958389365, 2.1583263863852635, 0.8168095027463367, 3.3557907081369653, 4.151603816416468, 2.7990053199341562 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] }, "processed": {"acEffectiveFrequency": 707644.9123862737, "average": -1.6666666666666643, "dutyCycle": 0.2421874999999996, "effectiveFrequency": 706482.814815305, "label": "Unipolar Rectangular", "offset": -53.33333333333333, "peak": 213.33333333333331, "peakToPeak": 266.66666666666663, "rms": 113.33333333333348, "thd": 0.9813747916705575 }, "waveform": {"ancillaryLabel": null, "data": [-53.33333333333333, 213.33333333333331, 213.33333333333331, -53.33333333333333, -53.33333333333333 ], "numberPeriods": null, "time": [0, 0, 0.000002421874999999996, 0.000002421874999999996, 0.00001 ] } } } ] } ] })");
        OpenMagnetics::MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 6);
        CHECK(masMagnetics.size() > 1);

        // std::cout << masMagnetics.size() << std::endl;
        // for (auto [mas, scoring] : masMagnetics) {
        //     std::string name = mas.get_magnetic().get_manufacturer_info().value().get_reference().value();
        //     auto masMagnetic = mas;
        //     OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
        //     auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
        //     auto outFile = outputFilePath;
        //     name = std::filesystem::path(std::regex_replace(std::string(name), std::regex(" "), "_")).string();
        //     name = std::filesystem::path(std::regex_replace(std::string(name), std::regex("\\,"), "_")).string();
        //     name = std::filesystem::path(std::regex_replace(std::string(name), std::regex("\\."), "_")).string();
        //     name = std::filesystem::path(std::regex_replace(std::string(name), std::regex("\\:"), "_")).string();
        //     name = std::filesystem::path(std::regex_replace(std::string(name), std::regex("\\/"), "_")).string();
        //     std::string filename = "MagneticAdviser_" + name + ".svg";
        //     std::cout << name << std::endl;
        //     outFile.append(filename);
        //     OpenMagnetics::Painter painter(outFile);

        //     painter.paint_core(masMagnetic.get_mutable_magnetic());
        //     painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
        //     painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
        //     painter.export_svg();
        // }
    }

    TEST(Test_MagneticAdviser_Web_2) {
        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->set_use_only_cores_in_stock(false);

        OpenMagnetics::InputsWrapper inputs = json::parse(R"({"designRequirements":{"magnetizingInductance":{"nominal":0.00001},"name":"My Design Requirements","turnsRatios":[]},"operatingPoints":[{"name":"Operating Point No. 1","conditions":{"ambientTemperature":42},"excitationsPerWinding":[{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.0000025,0.00001]},"processed":{"dutyCycle":0.25,"peakToPeak":10,"offset":0,"label":"Triangular","acEffectiveFrequency":128062.87364573497,"effectiveFrequency":128062.87364573497,"peak":4.999999999999998,"rms":2.887690890756071,"thd":0.37655984422983674},"harmonics":{"amplitudes":[1.057834375650657e-14,3.8218284740811654,1.3520347077861394,0.4253304394765621,2.0451658656097295e-15,0.15361202132455717,0.1511952491356625,0.07875264997066955,9.147290176737562e-16,0.04794839995573495,0.055136232876625904,0.03235759328251894,5.942500669304174e-16,0.023392967261331943,0.028681570420819563,0.017771008577657506,4.480960813987701e-16,0.014016396066471289,0.017807158306087034,0.011386475028370353,4.0120207228235354e-16,0.009474348163803114,0.012316262725232938,0.008041926341760212,3.291294404519727e-16,0.006942375707716075,0.00917327900077013,0.006081161484994139,1.7279005517465273e-16,0.005395200034140711,0.00721789399528569,0.00484110135817692,2.1149402582199242e-16,0.004388239086959991,0.005929248667932675,0.00401450896779898,2.157255682812684e-16,0.003703646967008042,0.005045721135369252,0.0034434847131706746,1.6787221123065324e-16,0.0032247860657256414,0.004424650744237714,0.0030404540187279378,2.3125579070802604e-16,0.0028849794948139283,0.003983384372505707,0.0027540516952747904,1.4216315610418345e-16,0.0026442785821521567,0.0036719561660328304,0.002552983887311765,1.889676895636462e-16,0.0024780582850305525,0.0034594522045110314,0.0024178496118518887,3.358503416551359e-16,0.002371081787894119,0.003326834522131644,0.0023367953000156096,1.7281675154410565e-16,0.0023143043115014237,0.0032630646073169567,0.002303167019442115],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}},"voltage":{"waveform":{"data":[-2.5,7.5,7.5,-2.5,-2.5],"time":[0,0,0.0000025,0.0000025,0.00001]},"processed":{"dutyCycle":0.25,"peakToPeak":10,"offset":0,"label":"Rectangular","acEffectiveFrequency":690362.0785055,"effectiveFrequency":690304.6703245763,"peak":7.5,"rms":4.284784125250653,"thd":0.9507084995254541},"harmonics":{"amplitudes":[0.078125,4.390192203074813,3.1805418164042503,1.6082999269932894,0.15625,0.7853068032814723,1.0533519383460939,0.747228683429235,0.15625,0.38152854838936,0.6237849662140766,0.5097284946091541,0.15625,0.22395874930175205,0.4366894957016387,0.3969081532881607,0.15625,0.1387657119725448,0.3303628683669758,0.32998229455569084,0.15625,0.08454863298834173,0.2606873758724238,0.2849610019773425,0.15625,0.046391828940503435,0.21067873648230084,0.25205943526842706,0.15625,0.017597888807606342,0.1723953087083574,0.22653191551027138,0.15625,0.005294556301253328,0.14161674515924202,0.20579082584016142,0.15625,0.02426175514926864,0.11588289785500591,0.18829805684603593,0.15625,0.040521304735682134,0.09365264588780098,0.17307466014248524,0.15625,0.05487174035875039,0.07390074623301954,0.1594602594300079,0.15625,0.06786649949748992,0.05590714395539467,0.1469848919856514,0.15625,0.07990999700987583,0.03913858752989208,0.13529576856110245,0.15625,0.09131439305868466,0.023177498052867507,0.1241125131778793,0.15625,0.10233553928853857,0.007676070276481761,0.11319769974196572],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}}}]}]})");
        OpenMagnetics::MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 6);
        CHECK(masMagnetics.size() > 1);

        // std::cout << masMagnetics.size() << std::endl;
        // for (auto [mas, scoring] : masMagnetics) {
        //     std::string name = mas.get_magnetic().get_manufacturer_info().value().get_reference().value();
        //     auto masMagnetic = mas;
        //     OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
        //     auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
        //     auto outFile = outputFilePath;
        //     name = std::filesystem::path(std::regex_replace(std::string(name), std::regex(" "), "_")).string();
        //     name = std::filesystem::path(std::regex_replace(std::string(name), std::regex("\\,"), "_")).string();
        //     name = std::filesystem::path(std::regex_replace(std::string(name), std::regex("\\."), "_")).string();
        //     name = std::filesystem::path(std::regex_replace(std::string(name), std::regex("\\:"), "_")).string();
        //     name = std::filesystem::path(std::regex_replace(std::string(name), std::regex("\\/"), "_")).string();
        //     std::string filename = "MagneticAdviser_" + name + ".svg";
        //     std::cout << name << std::endl;
        //     outFile.append(filename);
        //     OpenMagnetics::Painter painter(outFile);

        //     painter.paint_core(masMagnetic.get_mutable_magnetic());
        //     painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
        //     painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
        //     painter.export_svg();
        // }
    }

    TEST(Test_MagneticAdviser_Inductor) {
        srand (time(NULL));

        std::vector<double> turnsRatios;

        std::vector<OpenMagnetics::IsolationSide> isolationSides = {OpenMagnetics::IsolationSide::PRIMARY};

        double frequency = 50000;
        double peakToPeak = 90;
        double magnetizingInductance = 0.000110;
        double dutyCycle = 0.5;
        double dcCurrent = 0;
        double temperature = 25;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;

        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process_waveforms();

        OpenMagnetics::MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 1);

        if (masMagnetics.size() > 0) {
            auto masMagnetic = masMagnetics[0].first;
            OpenMagnetics::MagneticAdviser::preview_magnetic(masMagnetic);
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_MagneticAdviser_Inductor.svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagnetic.get_mutable_magnetic());
            painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            painter.export_svg();
        }
    }

    TEST(Test_MagneticAdviser_Inductor_Only_Toroidal_Cores) {
        srand (time(NULL));
        settings->reset();
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        settings->set_use_concentric_cores(false);


        std::vector<double> turnsRatios;

        std::vector<OpenMagnetics::IsolationSide> isolationSides = {OpenMagnetics::IsolationSide::PRIMARY};

        double frequency = 50000;
        double peakToPeak = 90;
        double magnetizingInductance = 0.000110;
        double dutyCycle = 0.5;
        double dcCurrent = 0;
        double temperature = 25;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;

        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process_waveforms();

        OpenMagnetics::MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 1);

        CHECK(masMagnetics.size() > 0);

        for (auto [masMagnetic, scoring] : masMagnetics) {
            OpenMagnetics::MagneticAdviser::preview_magnetic(masMagnetic);
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_MagneticAdviser_Inductor_Only_Toroidal_Cores.svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagnetic.get_mutable_magnetic());
            // painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            painter.export_svg();
        }
    }

    TEST(Test_MagneticAdviser_Inductor_PFC) {
        srand (time(NULL));
        settings->reset();
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        settings->set_use_concentric_cores(false);

        std::string file_path = __FILE__;
        auto inventory_path = file_path.substr(0, file_path.rfind("/")).append("/testData/pfc_current_waveform.ndjson");
        std::ifstream ndjsonFile(inventory_path);
        std::string jsonLine;


        std::vector<double> turnsRatios;

        std::vector<OpenMagnetics::IsolationSide> isolationSides = {OpenMagnetics::IsolationSide::PRIMARY,
                                                                    OpenMagnetics::IsolationSide::SECONDARY};
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

        // This chunk should go into a "import_csv" function in the future
        OpenMagnetics::Waveform waveform(windingExcitation["current"]["waveform"]);
        double frequency = 1.0 / (waveform.get_time().value().back() - waveform.get_time().value().front());
        auto sampledCurrentWaveform = OpenMagnetics::InputsWrapper::calculate_sampled_waveform(waveform, frequency);
        auto harmonics = OpenMagnetics::InputsWrapper::calculate_harmonics_data(sampledCurrentWaveform, frequency);
        settings->set_inputs_number_points_sampled_waveforms(2 * OpenMagnetics::round_up_size_to_power_of_2(harmonics.get_frequencies().back() / frequency));

        auto reconstructedWaveform = OpenMagnetics::InputsWrapper::reconstruct_signal(harmonics, frequency);
        settings->set_coil_allow_margin_tape(true);
        to_json(windingExcitation["current"]["waveform"], reconstructedWaveform);
        // ^^ This chunk

        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");

        auto outFile = outputFilePath;
        outFile.append("Test_MagneticAdviser_Inductor_PFC.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        painter.paint_waveform(reconstructedWaveform);
        // painter.paint_waveform(currentWaveform);
        painter.export_svg();

        operatingPoint["excitationsPerWinding"] = json::array();
        operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
        operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
        inputsJson["operatingPoints"].push_back(operatingPoint);

        inputsJson["designRequirements"] = json();
        inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 110e-6;
        inputsJson["designRequirements"]["turnsRatios"] = json::array();
        json turnsRatio;
        turnsRatio["nominal"] = 1;
        inputsJson["designRequirements"]["turnsRatios"].push_back(turnsRatio);
        OpenMagnetics::InputsWrapper inputs(inputsJson);

        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        inputs.process_waveforms();
        inputs.process_waveforms();
        OpenMagnetics::MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 3);

        CHECK(masMagnetics.size() > 0);

        for (auto [masMagnetic, scoring] : masMagnetics) {
            OpenMagnetics::MagneticAdviser::preview_magnetic(masMagnetic);
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_mutable_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_MagneticAdviser_Inductor_Only_Toroidal_Cores_" + std::to_string(scoring) + ".svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::SCATTER);
            settings->set_painter_logarithmic_scale(false);
            settings->set_painter_include_fringing(true);
            settings->set_painter_maximum_value_colorbar(std::nullopt);
            settings->set_painter_minimum_value_colorbar(std::nullopt);
            // settings->set_painter_number_points_x(200);
            // settings->set_painter_number_points_y(200);
            painter.paint_magnetic_field(masMagnetic.get_mutable_inputs().get_operating_point(0), masMagnetic.get_mutable_magnetic());
            painter.paint_core(masMagnetic.get_mutable_magnetic());
            // painter.paint_bobbin(masMagnetic.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetic.get_mutable_magnetic());
            painter.export_svg();
        }
    }

}