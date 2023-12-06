#include "Insulation.h"
#include <UnitTest++.h>
#include "TestingUtils.h"


SUITE(Insulation) {

    auto standardCoordinator = OpenMagnetics::InsulationCoordinator();
    auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
    auto cti = OpenMagnetics::Cti::GROUP_I;
    double maximumVoltageRms = 666;
    double maximumVoltagePeak = 800;
    double frequency = 30000;
    OpenMagnetics::DimensionWithTolerance altitude;
    OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;

    TEST(IEC_60664_Load_Data) {
        std::string dataString = R"({"IEC_60664-5": {"Table 2": {"inhomogeneusField": [[330, 0.00001], [400, 0.00002], [500, 0.00004], [600, 0.00006], [800, 0.0001], [1000, 0.00015], [1200, 0.00025], [1500, 0.0005], [2000, 0.001], [2500, 0.0015], [3000, 0.002], [4000, 0.0012], [5000, 0.0015], [6000, 0.002]], "homogeneusField": [[330, 0.00001], [400, 0.00002], [500, 0.00004], [600, 0.00006], [800, 0.0001], [1000, 0.00015], [1200, 0.0002], [1500, 0.0003], [2000, 0.00045], [2500, 0.0006], [3000, 0.0008]] }, "Table 3": {"inhomogeneusField": [[40, 0.000001], [60, 0.000002], [100, 0.000003], [120, 0.000004], [150, 0.000005], [200, 0.000006], [250, 0.000008], [330, 0.00001], [400, 0.00002], [500, 0.00004], [600, 0.00006], [800, 0.00013], [1000, 0.00026], [1200, 0.00042], [1500, 0.00076], [2000, 0.00127], [2500, 0.0018], [3000, 0.0024], [4000, 0.0012], [5000, 0.0015], [6000, 0.002]], "homogeneusField": [[40, 0.000001], [60, 0.000002], [100, 0.000003], [120, 0.000004], [150, 0.000005], [200, 0.000006], [250, 0.000008], [330, 0.00001], [400, 0.00002], [500, 0.00004], [600, 0.00006], [800, 0.0001], [1000, 0.00015], [1200, 0.0002], [1500, 0.0003], [2000, 0.00045], [2500, 0.0006], [3000, 0.0008]] }, "Table 4": {"P1": {"GROUP_I": [[40, 0.000025], [50, 0.000025], [63, 0.00004], [80, 0.000063], [100, 0.0001], [125, 0.00016], [160, 0.00025], [200, 0.0004], [250, 0.00056], [320, 0.00075], [400, 0.001], [500, 0.0013], [630, 0.0018], [800, 0.0024]], "GROUP_II": [[40, 0.000025], [50, 0.000025], [63, 0.00004], [80, 0.000063], [100, 0.0001], [125, 0.00016], [160, 0.00025], [200, 0.0004], [250, 0.00056], [320, 0.00075], [400, 0.001], [500, 0.0013], [630, 0.0018], [800, 0.0024]], "GROUP_IIIA": [[40, 0.000025], [50, 0.000025], [63, 0.00004], [80, 0.000063], [100, 0.0001], [125, 0.00016], [160, 0.00025], [200, 0.0004], [250, 0.00056], [320, 0.00075], [400, 0.001], [500, 0.0013], [630, 0.0018], [800, 0.0024]], "GROUP_IIIB": [[40, 0.000025], [50, 0.000025], [63, 0.00004], [80, 0.000063], [100, 0.0001], [125, 0.00016], [160, 0.00025], [200, 0.0004], [250, 0.00056], [320, 0.00075], [400, 0.001], [500, 0.0013], [630, 0.0018], [800, 0.0024]] }, "P2": {"GROUP_I": [[40, 0.00004], [50, 0.00004], [63, 0.000063], [80, 0.0001], [100, 0.00016], [125, 0.00025], [160, 0.0004], [200, 0.00063], [250, 0.001], [320, 0.0016], [400, 0.002]], "GROUP_II": [[40, 0.00004], [50, 0.00004], [63, 0.000063], [80, 0.0001], [100, 0.00016], [125, 0.00025], [160, 0.0004], [200, 0.00063], [250, 0.001], [320, 0.0016], [400, 0.002]], "GROUP_IIIA": [[40, 0.00004], [50, 0.00004], [63, 0.000063], [80, 0.0001], [100, 0.00016], [125, 0.00025], [160, 0.0004], [200, 0.00063], [250, 0.001], [320, 0.0016], [400, 0.002]] }, "P3": {"GROUP_I": [[40, 0.001], [50, 0.001], [63, 0.001], [80, 0.001], [100, 0.00125], [125, 0.0016], [160, 0.002]], "GROUP_II": [[40, 0.001], [50, 0.001], [63, 0.001], [80, 0.0011], [100, 0.0014], [125, 0.0018], [160, 0.0022]], "GROUP_IIIA": [[40, 0.001], [50, 0.001], [63, 0.001], [80, 0.00125], [100, 0.0016], [125, 0.002]] } } }, "IEC_60664-4": {"Table 1": [[600, 0.000065], [800, 0.00018], [1000, 0.0005], [1200, 0.0014], [1400, 0.00235], [1600, 0.004], [1800, 0.0067], [2000, 0.011]], "Table 2": {"100000": [[100, 0.000017], [200, 0.000042], [300, 0.000083], [400, 0.000125], [500, 0.000183], [600, 0.000267], [700, 0.000358], [800, 0.000450], [900, 0.000525], [1000, 0.000600], [1100, 0.000683], [1200, 0.000850], [1300, 0.001200], [1400, 0.001650], [1500, 0.002300], [1600, 0.003150], [1700, 0.004400], [1800, 0.006100]], "200000": [[300, 0.000090], [400, 0.000130], [500, 0.000190], [600, 0.000270], [700, 0.000380], [800, 0.000550], [900, 0.000820], [1000, 0.001150], [1100, 0.001700], [1200, 0.002400], [1300, 0.003500], [1400, 0.005000], [1500, 0.007300]], "400000": [[300, 0.000090], [400, 0.000150], [500, 0.000250], [600, 0.000400], [700, 0.000680], [800, 0.001100], [900, 0.001900], [1000, 0.003000], [1100, 0.005000], [1200, 0.008200]], "700000": [[300, 0.000090], [400, 0.000190], [500, 0.000400], [600, 0.000850], [700, 0.001900], [800, 0.003800], [900, 0.008700], [1000, 0.018000]], "1000000": [[300, 0.000090], [400, 0.000350], [500, 0.001500], [600, 0.005000], [700, 0.020000]], "2000000": [[200, 0.000150], [300, 0.000800], [400, 0.004500], [500, 0.020000]], "3000000": [[100, 0.000300], [200, 0.002800], [300, 0.020000]] } }, "IEC_60664-1": {"A.2": [[2000, 1.00], [3000, 1.14], [4000, 1.29], [5000, 1.48], [6000, 1.70], [7000, 1.95], [8000, 2.25], [9000, 2.62], [10000, 3.02], [15000, 6.67], [20000, 14.5]], "F.1": {"OVC_I": [[50, 330], [100, 500], [150, 800], [300, 1500], [600, 2500], [1000, 4000], [1250, 4000], [1500, 6000]], "OVC_II": [[50, 500], [100, 800], [150, 1500], [300, 2500], [600, 4000], [1000, 6000], [1250, 6000], [1500, 8000]], "OVC_III": [[50, 800], [100, 1500], [150, 2500], [300, 4000], [600, 6000], [1000, 8000], [1250, 8000], [1500, 10000]], "OVC_IV": [[50, 1500], [100, 2500], [150, 4000], [300, 6000], [600, 8000], [1000, 12000], [1250, 12000], [1500, 15000]] }, "F.2": {"inhomogeneusField": {"P1": [[330, 0.00001], [400, 0.00002], [500, 0.00004], [600, 0.00006], [800, 0.0001], [1000, 0.00015], [1200, 0.00025], [1500, 0.0005], [2000, 0.001], [2500, 0.0015], [3000, 0.002], [4000, 0.003], [5000, 0.004], [6000, 0.0055], [8000, 0.008], [10000, 0.011], [12000, 0.014], [15000, 0.018], [20000, 0.025], [25000, 0.033], [30000, 0.04], [40000, 0.06], [50000, 0.075], [60000, 0.09], [80000, 0.13], [100000, 0.17]], "P2": [[330, 0.0002], [400, 0.0002], [500, 0.0002], [600, 0.0002], [800, 0.0002], [1000, 0.0002], [1200, 0.00025], [1500, 0.0005], [2000, 0.001], [2500, 0.0015], [3000, 0.002], [4000, 0.003], [5000, 0.004], [6000, 0.0055], [8000, 0.008], [10000, 0.011], [12000, 0.014], [15000, 0.018], [20000, 0.025], [25000, 0.033], [30000, 0.04], [40000, 0.06], [50000, 0.075], [60000, 0.09], [80000, 0.13], [100000, 0.17]], "P3": [[330, 0.0008], [400, 0.0008], [500, 0.0008], [600, 0.0008], [800, 0.0008], [1000, 0.0008], [1200, 0.0008], [1500, 0.0008], [2000, 0.001], [2500, 0.0015], [3000, 0.002], [4000, 0.003], [5000, 0.004], [6000, 0.0055], [8000, 0.008], [10000, 0.011], [12000, 0.014], [15000, 0.018], [20000, 0.025], [25000, 0.033], [30000, 0.04], [40000, 0.06], [50000, 0.075], [60000, 0.09], [80000, 0.13], [100000, 0.17]] }, "homogeneusField": {"P1": [[330, 0.00001], [400, 0.00002], [500, 0.00004], [600, 0.00006], [800, 0.0001], [1000, 0.00015], [1200, 0.0002], [1500, 0.0003], [2000, 0.00045], [2500, 0.0006], [3000, 0.0008], [4000, 0.0012], [5000, 0.0015], [6000, 0.002], [8000, 0.003], [10000, 0.0035], [12000, 0.0045], [15000, 0.0055], [20000, 0.008], [25000, 0.01], [30000, 0.0125], [40000, 0.017], [50000, 0.022], [60000, 0.027], [80000, 0.035], [100000, 0.045]], "P2": [[330, 0.0002], [400, 0.0002], [500, 0.0002], [600, 0.0002], [800, 0.0002], [1000, 0.0002], [1200, 0.0002], [1500, 0.0003], [2000, 0.00045], [2500, 0.0006], [3000, 0.0008], [4000, 0.0012], [5000, 0.0015], [6000, 0.002], [8000, 0.003], [10000, 0.0035], [12000, 0.0045], [15000, 0.0055], [20000, 0.008], [25000, 0.01], [30000, 0.0125], [40000, 0.017], [50000, 0.022], [60000, 0.027], [80000, 0.035], [100000, 0.045]], "P3": [[330, 0.0008], [400, 0.0008], [500, 0.0008], [600, 0.0008], [800, 0.0008], [1000, 0.0008], [1200, 0.0008], [1500, 0.0008], [2000, 0.0008], [2500, 0.0008], [3000, 0.0008], [4000, 0.0012], [5000, 0.0015], [6000, 0.002], [8000, 0.003], [10000, 0.0035], [12000, 0.0045], [15000, 0.0055], [20000, 0.008], [25000, 0.01], [30000, 0.0125], [40000, 0.017], [50000, 0.022], [60000, 0.027], [80000, 0.035], [100000, 0.045]] } }, "F.3": [[12.5, 12.5], [25, 25], [30, 32], [42, 50], [48, 50], [50, 50], [60, 63], [100, 100], [110, 125], [120, 125], [150, 160], [200, 200], [220, 250], [300, 320], [600, 630], [1000, 1000], [1500, 1500]], "F.5": {"PRINTED": {"P1": {"GROUP_I": [[10, 0.000025], [12.5, 0.000025], [16, 0.000025], [20, 0.000025], [25, 0.000025], [32, 0.000025], [40, 0.000025], [50, 0.000025], [63, 0.000040], [80, 0.000063], [100, 0.000100], [125, 0.000160], [160, 0.000250], [200, 0.000400], [250, 0.000560], [320, 0.000750], [400, 0.001000], [500, 0.001300], [630, 0.001800], [800, 0.002400], [1000, 0.003200]], "GROUP_II": [[10, 0.000025], [12.5, 0.000025], [16, 0.000025], [20, 0.000025], [25, 0.000025], [32, 0.000025], [40, 0.000025], [50, 0.000025], [63, 0.000040], [80, 0.000063], [100, 0.000100], [125, 0.000160], [160, 0.000250], [200, 0.000400], [250, 0.000560], [320, 0.000750], [400, 0.001000], [500, 0.001300], [630, 0.001800], [800, 0.002400], [1000, 0.003200]], "GROUP_IIIA": [[10, 0.000025], [12.5, 0.000025], [16, 0.000025], [20, 0.000025], [25, 0.000025], [32, 0.000025], [40, 0.000025], [50, 0.000025], [63, 0.000040], [80, 0.000063], [100, 0.000100], [125, 0.000160], [160, 0.000250], [200, 0.000400], [250, 0.000560], [320, 0.000750], [400, 0.001000], [500, 0.001300], [630, 0.001800], [800, 0.002400], [1000, 0.003200]], "GROUP_IIIB": [[10, 0.000025], [12.5, 0.000025], [16, 0.000025], [20, 0.000025], [25, 0.000025], [32, 0.000025], [40, 0.000025], [50, 0.000025], [63, 0.000040], [80, 0.000063], [100, 0.000100], [125, 0.000160], [160, 0.000250], [200, 0.000400], [250, 0.000560], [320, 0.000750], [400, 0.001000], [500, 0.001300], [630, 0.001800], [800, 0.002400], [1000, 0.003200]] }, "P2": {"GROUP_I": [[10, 0.000040], [12.5, 0.000040], [16, 0.000040], [20, 0.000040], [25, 0.000040], [32, 0.000040], [40, 0.000040], [50, 0.000040], [63, 0.000063], [80, 0.000100], [100, 0.000160], [125, 0.000250], [160, 0.000400], [200, 0.000630], [250, 0.001000], [320, 0.001600], [400, 0.002000], [500, 0.002500], [630, 0.003200], [800, 0.004000], [1000, 0.005000]], "GROUP_II": [[10, 0.000040], [12.5, 0.000040], [16, 0.000040], [20, 0.000040], [25, 0.000040], [32, 0.000040], [40, 0.000040], [50, 0.000040], [63, 0.000063], [80, 0.000100], [100, 0.000160], [125, 0.000250], [160, 0.000400], [200, 0.000630], [250, 0.001000], [320, 0.001600], [400, 0.002000], [500, 0.002500], [630, 0.003200], [800, 0.004000], [1000, 0.005000]], "GROUP_IIIA": [[10, 0.000040], [12.5, 0.000040], [16, 0.000040], [20, 0.000040], [25, 0.000040], [32, 0.000040], [40, 0.000040], [50, 0.000040], [63, 0.000063], [80, 0.000100], [100, 0.000160], [125, 0.000250], [160, 0.000400], [200, 0.000630], [250, 0.001000], [320, 0.001600], [400, 0.002000], [500, 0.002500], [630, 0.003200], [800, 0.004000], [1000, 0.005000]] } }, "WOUND": {"P1": {"GROUP_I": [[10, 0.000080], [12.5, 0.000090], [16, 0.000100], [20, 0.000110], [25, 0.000125], [32, 0.000140], [40, 0.000160], [50, 0.000180], [63, 0.000200], [80, 0.000220], [100, 0.000250], [125, 0.000280], [160, 0.000320], [200, 0.000420], [250, 0.000560], [320, 0.000750], [400, 0.001000], [500, 0.001300], [630, 0.001800], [800, 0.002400], [1000, 0.003200], [1250, 0.004200], [1600, 0.005600], [2000, 0.007500], [2500, 0.010000], [3200, 0.012500], [4000, 0.016], [5000, 0.02], [6300, 0.025], [8000, 0.032], [10000, 0.04], [12500, 0.05], [16000, 0.063], [20000, 0.08], [25000, 0.1], [32000, 0.125], [40000, 0.16], [50000, 0.2], [63000, 0.25]], "GROUP_II": [[10, 0.000080], [12.5, 0.000090], [16, 0.000100], [20, 0.000110], [25, 0.000125], [32, 0.000140], [40, 0.000160], [50, 0.000180], [63, 0.000200], [80, 0.000220], [100, 0.000250], [125, 0.000280], [160, 0.000320], [200, 0.000420], [250, 0.000560], [320, 0.000750], [400, 0.001000], [500, 0.001300], [630, 0.001800], [800, 0.002400], [1000, 0.003200], [1250, 0.004200], [1600, 0.005600], [2000, 0.007500], [2500, 0.010000], [3200, 0.012500], [4000, 0.016], [5000, 0.02], [6300, 0.025], [8000, 0.032], [10000, 0.04], [12500, 0.05], [16000, 0.063], [20000, 0.08], [25000, 0.1], [32000, 0.125], [40000, 0.16], [50000, 0.2], [63000, 0.25]], "GROUP_IIIA": [[10, 0.000080], [12.5, 0.000090], [16, 0.000100], [20, 0.000110], [25, 0.000125], [32, 0.000140], [40, 0.000160], [50, 0.000180], [63, 0.000200], [80, 0.000220], [100, 0.000250], [125, 0.000280], [160, 0.000320], [200, 0.000420], [250, 0.000560], [320, 0.000750], [400, 0.001000], [500, 0.001300], [630, 0.001800], [800, 0.002400], [1000, 0.003200], [1250, 0.004200], [1600, 0.005600], [2000, 0.007500], [2500, 0.010000], [3200, 0.012500], [4000, 0.016], [5000, 0.02], [6300, 0.025], [8000, 0.032], [10000, 0.04], [12500, 0.05], [16000, 0.063], [20000, 0.08], [25000, 0.1], [32000, 0.125], [40000, 0.16], [50000, 0.2], [63000, 0.25]], "GROUP_IIIB": [[10, 0.000080], [12.5, 0.000090], [16, 0.000100], [20, 0.000110], [25, 0.000125], [32, 0.000140], [40, 0.000160], [50, 0.000180], [63, 0.000200], [80, 0.000220], [100, 0.000250], [125, 0.000280], [160, 0.000320], [200, 0.000420], [250, 0.000560], [320, 0.000750], [400, 0.001000], [500, 0.001300], [630, 0.001800], [800, 0.002400], [1000, 0.003200], [1250, 0.004200], [1600, 0.005600], [2000, 0.007500], [2500, 0.010000], [3200, 0.012500], [4000, 0.016], [5000, 0.02], [6300, 0.025], [8000, 0.032], [10000, 0.04], [12500, 0.05], [16000, 0.063], [20000, 0.08], [25000, 0.1], [32000, 0.125], [40000, 0.16], [50000, 0.2], [63000, 0.25]] }, "P2": {"GROUP_I": [[10, 0.00040], [12.5, 0.00042], [16, 0.00045], [20, 0.00048], [25, 0.00050], [32, 0.00053], [40, 0.00056], [50, 0.00060], [63, 0.00063], [80, 0.00067], [100, 0.00071], [125, 0.00075], [160, 0.00080], [200, 0.00100], [250, 0.00125], [320, 0.00160], [400, 0.00200], [500, 0.00250], [630, 0.00320], [800, 0.00400], [1000, 0.00500], [1250, 0.00630], [1600, 0.00800], [2000, 0.01000], [2500, 0.01250], [3200, 0.01600], [4000, 0.02], [5000, 0.025], [6300, 0.032], [8000, 0.04], [10000, 0.05], [12500, 0.063], [16000, 0.08], [20000, 0.1], [25000, 0.125], [32000, 0.16], [40000, 0.2], [50000, 0.25], [63000, 0.32]], "GROUP_II": [[10, 0.000400], [12.5, 0.000420], [16, 0.000450], [20, 0.000480], [25, 0.000500], [32, 0.000530], [40, 0.000800], [50, 0.000850], [63, 0.000900], [80, 0.000950], [100, 0.001000], [125, 0.001050], [160, 0.001100], [200, 0.001400], [250, 0.001800], [320, 0.002200], [400, 0.002800], [500, 0.003600], [630, 0.004500], [800, 0.005600], [1000, 0.007100], [1250, 0.009000], [1600, 0.011000], [2000, 0.014000], [2500, 0.018000], [3200, 0.022000], [4000, 0.028], [5000, 0.036], [6300, 0.045], [8000, 0.056], [10000, 0.071], [12500, 0.09], [16000, 0.11], [20000, 0.14], [25000, 0.18], [32000, 0.22], [40000, 0.28], [50000, 0.36], [63000, 0.45]], "GROUP_IIIA": [[10, 0.00040], [12.5, 0.00042], [16, 0.00045], [20, 0.00048], [25, 0.00050], [32, 0.00053], [40, 0.00110], [50, 0.00120], [63, 0.00125], [80, 0.00130], [100, 0.00140], [125, 0.00150], [160, 0.00160], [200, 0.00200], [250, 0.00250], [320, 0.00320], [400, 0.00400], [500, 0.00500], [630, 0.00630], [800, 0.00800], [1000, 0.01000], [1250, 0.01250], [1600, 0.01600], [2000, 0.02000], [2500, 0.02500], [3200, 0.03200], [4000, 0.04], [5000, 0.05], [6300, 0.063], [8000, 0.08], [10000, 0.1], [12500, 0.125], [16000, 0.16], [20000, 0.2], [25000, 0.25], [32000, 0.32], [40000, 0.4], [50000, 0.5], [63000, 0.6]], "GROUP_IIIB": [[10, 0.00040], [12.5, 0.00042], [16, 0.00045], [20, 0.00048], [25, 0.00050], [32, 0.00053], [40, 0.00110], [50, 0.00120], [63, 0.00125], [80, 0.00130], [100, 0.00140], [125, 0.00150], [160, 0.00160], [200, 0.00200], [250, 0.00250], [320, 0.00320], [400, 0.00400], [500, 0.00500], [630, 0.00630], [800, 0.00800], [1000, 0.01000], [1250, 0.01250], [1600, 0.01600], [2000, 0.02000], [2500, 0.02500], [3200, 0.03200], [4000, 0.04], [5000, 0.05], [6300, 0.063], [8000, 0.08], [10000, 0.1], [12500, 0.125], [16000, 0.16], [20000, 0.2], [25000, 0.25], [32000, 0.32], [40000, 0.4], [50000, 0.5], [63000, 0.6]] }, "P3": {"GROUP_I": [[10, 0.00100], [12.5, 0.00105], [16, 0.00110], [20, 0.00120], [25, 0.00125], [32, 0.00130], [40, 0.00140], [50, 0.00150], [63, 0.00160], [80, 0.00170], [100, 0.00180], [125, 0.00190], [160, 0.00200], [200, 0.00250], [250, 0.00320], [320, 0.00400], [400, 0.00500], [500, 0.00630], [630, 0.00800], [800, 0.01000], [1000, 0.01250], [1250, 0.01600], [1600, 0.02000], [2000, 0.02500], [2500, 0.03200], [3200, 0.04000], [4000, 0.05], [5000, 0.063], [6300, 0.08], [8000, 0.1], [10000, 0.125]], "GROUP_II": [[10, 0.00100], [12.5, 0.00105], [16, 0.00110], [20, 0.00120], [25, 0.00125], [32, 0.00130], [40, 0.00160], [50, 0.00170], [63, 0.00180], [80, 0.00190], [100, 0.00200], [125, 0.00210], [160, 0.00220], [200, 0.00280], [250, 0.00360], [320, 0.00450], [400, 0.00560], [500, 0.00710], [630, 0.00900], [800, 0.01100], [1000, 0.01400], [1250, 0.01800], [1600, 0.02200], [2000, 0.02800], [2500, 0.03600], [3200, 0.04500], [4000, 0.056], [5000, 0.071], [6300, 0.09 ], [8000, 0.11 ], [10000, 0.14]], "GROUP_IIIA": [[10, 0.00100], [12.5, 0.00105], [16, 0.00110], [20, 0.00120], [25, 0.00125], [32, 0.00130], [40, 0.00180], [50, 0.00190], [63, 0.00200], [80, 0.00210], [100, 0.00220], [125, 0.00240], [160, 0.00250], [200, 0.00320], [250, 0.00400], [320, 0.00500], [400, 0.00630], [500, 0.00800], [630, 0.01000], [800, 0.01250], [1000, 0.01600], [1250, 0.02000], [1600, 0.02500], [2000, 0.03200], [2500, 0.04000], [3200, 0.05000], [4000, 0.063], [5000, 0.08], [6300, 0.1], [8000, 0.125], [10000, 0.16]], "GROUP_IIIB": [[10, 0.00100], [12.5, 0.00105], [16, 0.00110], [20, 0.00120], [25, 0.00125], [32, 0.00130], [40, 0.00180], [50, 0.00190], [63, 0.00200], [80, 0.00210], [100, 0.00220], [125, 0.00240], [160, 0.00250], [200, 0.00320], [250, 0.00400], [320, 0.00500], [400, 0.00630], [500, 0.00800], [630, 0.01000], [800, 0.01250], [1000, 0.01600], [1250, 0.02000], [1600, 0.02500], [2000, 0.03200], [2500, 0.04000], [3200, 0.05000], [4000, 0.063], [5000, 0.08], [6300, 0.1], [8000, 0.125], [10000, 0.16]] } } }, "F.8": {"inhomogeneusField": [[40, 0.000001], [60, 0.000002], [100, 0.000003], [120, 0.000004], [150, 0.000005], [200, 0.000006], [250, 0.000008], [330, 0.000010], [400, 0.000020], [500, 0.000040], [600, 0.000060], [800, 0.000130], [1000, 0.000260], [1200, 0.000420], [1500, 0.000760], [2000, 0.001270], [2500, 0.001800], [3000, 0.002400], [4000, 0.003800], [5000, 0.005700], [6000, 0.007900], [8000, 0.011000], [10000, 0.015200], [12000, 0.019000], [15000, 0.025000], [20000, 0.034000], [25000, 0.044000], [30000, 0.055000], [40000, 0.077000], [50000, 0.100000]], "homogeneusField": [[40, 0.000001], [60, 0.000002], [100, 0.000003], [120, 0.000004], [150, 0.000005], [200, 0.000006], [250, 0.000008], [330, 0.000010], [400, 0.000020], [500, 0.000040], [600, 0.000060], [800, 0.000100], [1000, 0.000150], [1200, 0.000200], [1500, 0.000300], [2000, 0.000450], [2500, 0.000600], [3000, 0.000800], [4000, 0.001200], [5000, 0.001500], [6000, 0.002000], [8000, 0.003000], [10000, 0.003500], [12000, 0.004500], [15000, 0.005500], [20000, 0.008000], [25000, 0.010000], [30000, 0.012500], [40000, 0.017000], [50000, 0.022000], [60000, 0.027000], [80000, 0.035000], [100000, 0.045000]] } } })";
        json data = json::parse(dataString);
        auto standard = OpenMagnetics::InsulationIEC60664Model(data);
        OpenMagnetics::DimensionWithTolerance altitude;
        altitude.set_maximum(2000);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
        double maximumVoltageRms = 666;
        double maximumVoltagePeak = 800;
        double frequency = 30000;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND); auto solidInsulation = standard.calculate_solid_insulation(inputs);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0024, creepageDistance);
    }

    TEST(Test_Coordinated_Creepage_Distance) {

        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641, OpenMagnetics::InsulationStandards::IEC_623681};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standardCoordinator.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0024, creepageDistance);
    }

    TEST(Test_Coordinated_Clearance) {
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641, OpenMagnetics::InsulationStandards::IEC_623681};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standardCoordinator.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }


}

SUITE(CreepageDistance_IEC_60664) {

    auto standard = OpenMagnetics::InsulationIEC60664Model();
    auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
    auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
    double maximumVoltageRms = 666;
    double maximumVoltagePeak = 800;
    double frequency = 30000;
    OpenMagnetics::DimensionWithTolerance altitude;
    OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;


    TEST(Creepage_Distance_Basic_P1_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0024, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0048, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.004, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.008, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.01, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.02, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0024, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0048, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0056, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0112, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.011, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.022, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0024, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0048, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.008, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.016, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0125, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.025, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_I_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK(0.0024 < creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_I_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK(0.0048 < creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_I_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK(0.004 < creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_I_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK(0.008 < creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_I_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK(0.01 < creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_I_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK(0.02 < creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_II_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK(0.0024 < creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_II_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK(0.0048 < creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_II_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK(0.0056 < creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_II_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK(0.0112 < creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_II_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK(0.011 < creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_II_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK(0.022 < creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_IIIA_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK(0.0024 < creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_IIIA_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK(0.0048 < creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_IIIA_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK(0.008 < creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_IIIA_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK(0.016 < creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_IIIA_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK(0.0125 < creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_IIIA_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK(0.025 < creepageDistance);
    }
}

SUITE(Clearance_IEC_60664) {
    auto standard = OpenMagnetics::InsulationIEC60664Model();
    auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
    auto cti = OpenMagnetics::Cti::GROUP_I;
    double maximumVoltageRms = 69;
    double maximumVoltagePeak = 260;
    double frequency = 30000;
    OpenMagnetics::DimensionWithTolerance altitude;
    OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;

    TEST(Clearance_Basic_P1_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.00004, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.00006, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0002, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0002, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0008, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0008, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0001, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.00015, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0002, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0002, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0008, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0008, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0005, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.001, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0005, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.001, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0008, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.001, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0015, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.002, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0015, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.002, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0015, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.002, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_I_High_Altitude_Low_Frequency) {
        altitude.set_maximum(8000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.00004 * 2.25, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_I_Low_Altitude_Low_Frequency_High_Voltage) {
        maximumVoltageRms = 666;
        maximumVoltagePeak = 800;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_I_Low_Altitude_High_Frequency_High_Voltage) {
        frequency = 500000;
        maximumVoltageRms = 666;
        maximumVoltagePeak = 800;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }
}

SUITE(CreepageDistance_IEC_62368) {

    auto standard = OpenMagnetics::InsulationIEC62368Model();
    auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_623681};
    auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
    double maximumVoltageRms = 666;
    double maximumVoltagePeak = 800;
    double frequency = 30000;
    OpenMagnetics::DimensionWithTolerance altitude;
    OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;


    TEST(Creepage_Distance_Basic_P1_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.002, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0039, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0034, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0068, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0085, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0169, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.002, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0039, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0048, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0095, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0095, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0189, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.002, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0039, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0067, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0134, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0106, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0211, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_I_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.002, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_I_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0039, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_I_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0034, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_I_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0068, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_I_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0085, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_I_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0169, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_II_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.002, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_II_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0039, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_II_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0048, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_II_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0095, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_II_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0095, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_II_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0189, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_IIIA_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.002, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_IIIA_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0039, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_IIIA_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0067, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_IIIA_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0134, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_IIIA_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0106, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_IIIA_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0211, creepageDistance);
    }
}

SUITE(Clearance_IEC_62368) {
    auto standard = OpenMagnetics::InsulationIEC62368Model();
    auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_623681};
        auto cti = OpenMagnetics::Cti::GROUP_I;
    double maximumVoltageRms = 666;
    double maximumVoltagePeak = 800;
    double frequency = 30000;
    OpenMagnetics::DimensionWithTolerance altitude;
    OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;

    TEST(Clearance_Basic_P1_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0018, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0036, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0018, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0036, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0018, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0036, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0018, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0036, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0018, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0036, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0018, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0036, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0055, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0055, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0055, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0055, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.008, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0055, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.008, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0055, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.008, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_I_High_Altitude_Low_Frequency) {
        altitude.set_maximum(5000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.00267, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_I_High_Altitude_Low_Frequency) {
        altitude.set_maximum(5000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.00533, clearance);
    }    

    TEST(Clearance_Basic_P1_OVC_I_Low_Altitude_High_Frequency_High_Voltage_Peak) {
        altitude.set_maximum(2000);
        frequency = 400000;
        maximumVoltagePeak = 2000;
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0106, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_I_Low_Altitude_High_Frequency_High_Voltage_Peak) {
        altitude.set_maximum(2000);
        frequency = 400000;
        maximumVoltagePeak = 2000;
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0212, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_I_Low_Altitude_High_Frequency_High_Voltage_Peak) {
        altitude.set_maximum(2000);
        frequency = 400000;
        maximumVoltagePeak = 2000;
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0132, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_I_Low_Altitude_High_Frequency_High_Voltage_Peak) {
        altitude.set_maximum(2000);
        frequency = 400000;
        maximumVoltagePeak = 2000;
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0264, clearance);
    }

    TEST(Clearance_Printed_Basic) {
        altitude.set_maximum(2000);
        frequency = 100000;
        maximumVoltagePeak = 2000;
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::PRINTED);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0038, clearance);
    }

    TEST(Clearance_Printed_Reinforced) {
        altitude.set_maximum(2000);
        frequency = 100000;
        maximumVoltagePeak = 2000;
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::PRINTED);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0044, clearance);
    }

}
