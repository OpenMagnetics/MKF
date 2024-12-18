#include "Constants.h"
#include "CoreLosses.h"
#include "Settings.h"
#include "InputsWrapper.h"
#include "MagnetizingInductance.h"
#include "Reluctance.h"
#include "TestingUtils.h"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

std::map<std::string, std::map<int, std::map<std::string, double>>> dynamicCoefficients = {
    {
        "3F4",
        {
            {1000000,
             {
                 {"alpha", 2.319945632800448},
                 {"beta", 2.326152651677128},
                 {"ct0", 0.000317296062135492},
                 {"ct1", 3.763091992620556e-06},
                 {"ct2", 3.521451072997681e-08},
                 {"k", 0.047520631299750656},
             }},
            {500000,
             {
                 {"alpha", 2.2270388480289673},
                 {"beta", 2.3370905115783476},
                 {"ct0", 0.0005483667014872101},
                 {"ct1", 6.822239047393187e-06},
                 {"ct2", 6.308173416993905e-08},
                 {"k", 0.09948590338409125},
             }},
            {800000,
             {
                 {"alpha", 2.319945632800448},
                 {"beta", 2.326152651677128},
                 {"ct0", 0.000317296062135492},
                 {"ct1", 3.763091992620556e-06},
                 {"ct2", 3.521451072997681e-08},
                 {"k", 0.047520631299750656},
             }},
            {300000,
             {
                 {"alpha", 1.8864984600467845},
                 {"beta", 2.6271189944168762},
                 {"ct0", 0.002350040699111803},
                 {"ct1", 2.9331101755357237e-05},
                 {"ct2", 2.2449305959577767e-07},
                 {"k", 8.252492522049222},
             }},
            {400000,
             {
                 {"alpha", 1.3817827522244912},
                 {"beta", 2.3596181724302174},
                 {"ct0", 0.0030695032438131215},
                 {"ct1", 4.062568887844067e-05},
                 {"ct2", 3.561155879040697e-07},
                 {"k", 1331.770034940966},
             }},
            {50000,
             {
                 {"alpha", 1.0393937536437856},
                 {"beta", 2.916512195344141},
                 {"ct0", 0.004600052615555796},
                 {"ct1", 5.3225563373133944e-05},
                 {"ct2", 3.557608906354458e-07},
                 {"k", 369318.91764390975},
             }},
            {25000,
             {
                 {"alpha", 1.0393937536437856},
                 {"beta", 2.916512195344141},
                 {"ct0", 0.004600052615555796},
                 {"ct1", 5.3225563373133944e-05},
                 {"ct2", 3.557608906354458e-07},
                 {"k", 369318.91764390975},
             }},
            {100000,
             {
                 {"alpha", 1.29372652797385},
                 {"beta", 2.9992082547433725},
                 {"ct0", 0.004201708719962025},
                 {"ct1", 5.195425637860509e-05},
                 {"ct2", 3.5001294769795715e-07},
                 {"k", 26361.71634950812},
             }},
            {200000,
             {
                 {"alpha", 1.1659757620484792},
                 {"beta", 2.809880734072415},
                 {"ct0", 0.004523095913359234},
                 {"ct1", 5.809673530388853e-05},
                 {"ct2", 4.0819693608071845e-07},
                 {"k", 64609.88337459773},
             }},
        },
    },
    {
        "N49",
        {
            {1000000,
             {
                 {"alpha", 1.2584328161249612},
                 {"beta", 2.443194133466379},
                 {"ct0", 0.004008571051875008},
                 {"ct1", 4.34285989724493e-05},
                 {"ct2", 5.56224171234087e-07},
                 {"k", 3859.6637037342953},
             }},
            {500000,
             {
                 {"alpha", 2.2476484059505495},
                 {"beta", 2.490359501523948},
                 {"ct0", 0.0029444317262005484},
                 {"ct1", 2.345183702579143e-05},
                 {"ct2", 3.387858472797546e-07},
                 {"k", 0.012595006528621926},
             }},
            {800000,
             {
                 {"alpha", 1.2584328161249612},
                 {"beta", 2.443194133466379},
                 {"ct0", 0.004008571051875008},
                 {"ct1", 4.34285989724493e-05},
                 {"ct2", 5.56224171234087e-07},
                 {"k", 3859.6637037342953},
             }},
            {300000,
             {
                 {"alpha", 1.6199173550878625},
                 {"beta", 2.906402440536201},
                 {"ct0", 0.0031619656170100794},
                 {"ct1", 4.1112538084232765e-05},
                 {"ct2", 4.942169715197847e-07},
                 {"k", 148.30054599084795},
             }},
            {400000,
             {
                 {"alpha", 3.104644728196987},
                 {"beta", 2.5919391447870597},
                 {"ct0", 1.859009254165758e-05},
                 {"ct1", 1.7087201755673073e-07},
                 {"ct2", 2.0528734308554635e-09},
                 {"k", 4.324316353158944e-05},
             }},
            {50000,
             {
                 {"alpha", 1.0034329035439096},
                 {"beta", 2.9538857892019634},
                 {"ct0", 0.0034755214933370124},
                 {"ct1", 5.1878796028138934e-05},
                 {"ct2", 5.431559753579025e-07},
                 {"k", 319245.2578685631},
             }},
            {25000,
             {
                 {"alpha", 1.0034329035439096},
                 {"beta", 2.9538857892019634},
                 {"ct0", 0.0034755214933370124},
                 {"ct1", 5.1878796028138934e-05},
                 {"ct2", 5.431559753579025e-07},
                 {"k", 319245.2578685631},
             }},
            {100000,
             {
                 {"alpha", 1.3345013134805503},
                 {"beta", 3.091014997668596},
                 {"ct0", 0.003074480332680605},
                 {"ct1", 5.314117675943326e-05},
                 {"ct2", 5.48645765579563e-07},
                 {"k", 10768.644284113128},
             }},
            {200000,
             {
                 {"alpha", 1.5352695641272587},
                 {"beta", 3.0268823360088155},
                 {"ct0", 0.002574021633028376},
                 {"ct1", 4.0405038232170986e-05},
                 {"ct2", 4.311439514950848e-07},
                 {"k", 842.5250133813021},
             }},
        },
    },
    {
        "N49",
        {
            {1000000,
             {
                 {"alpha", 1.2584328161249612},
                 {"beta", 2.443194133466379},
                 {"ct0", 0.004008571051875008},
                 {"ct1", 4.34285989724493e-05},
                 {"ct2", 5.56224171234087e-07},
                 {"k", 3859.6637037342953},
             }},
            {500000,
             {
                 {"alpha", 2.2476484059505495},
                 {"beta", 2.490359501523948},
                 {"ct0", 0.0029444317262005484},
                 {"ct1", 2.345183702579143e-05},
                 {"ct2", 3.387858472797546e-07},
                 {"k", 0.012595006528621926},
             }},
            {800000,
             {
                 {"alpha", 1.2584328161249612},
                 {"beta", 2.443194133466379},
                 {"ct0", 0.004008571051875008},
                 {"ct1", 4.34285989724493e-05},
                 {"ct2", 5.56224171234087e-07},
                 {"k", 3859.6637037342953},
             }},
            {300000,
             {
                 {"alpha", 1.6199173550878625},
                 {"beta", 2.906402440536201},
                 {"ct0", 0.0031619656170100794},
                 {"ct1", 4.1112538084232765e-05},
                 {"ct2", 4.942169715197847e-07},
                 {"k", 148.30054599084795},
             }},
            {400000,
             {
                 {"alpha", 3.104644728196987},
                 {"beta", 2.5919391447870597},
                 {"ct0", 1.859009254165758e-05},
                 {"ct1", 1.7087201755673073e-07},
                 {"ct2", 2.0528734308554635e-09},
                 {"k", 4.324316353158944e-05},
             }},
            {50000,
             {
                 {"alpha", 1.0034329035439096},
                 {"beta", 2.9538857892019634},
                 {"ct0", 0.0034755214933370124},
                 {"ct1", 5.1878796028138934e-05},
                 {"ct2", 5.431559753579025e-07},
                 {"k", 319245.2578685631},
             }},
            {25000,
             {
                 {"alpha", 1.0034329035439096},
                 {"beta", 2.9538857892019634},
                 {"ct0", 0.0034755214933370124},
                 {"ct1", 5.1878796028138934e-05},
                 {"ct2", 5.431559753579025e-07},
                 {"k", 319245.2578685631},
             }},
            {100000,
             {
                 {"alpha", 1.3345013134805503},
                 {"beta", 3.091014997668596},
                 {"ct0", 0.003074480332680605},
                 {"ct1", 5.314117675943326e-05},
                 {"ct2", 5.48645765579563e-07},
                 {"k", 10768.644284113128},
             }},
            {200000,
             {
                 {"alpha", 1.5352695641272587},
                 {"beta", 3.0268823360088155},
                 {"ct0", 0.002574021633028376},
                 {"ct1", 4.0405038232170986e-05},
                 {"ct2", 4.311439514950848e-07},
                 {"k", 842.5250133813021},
             }},
        },
    },
    {
        "3C94",
        {
            {1000000,
             {
                 {"alpha", 4.123853571168164},
                 {"beta", 2.07914339855886},
                 {"ct0", 6.052001156627106e-09},
                 {"ct1", 1.3731037996663199e-11},
                 {"ct2", 3.315293060549009e-14},
                 {"k", 9.812805086724211e-08},
             }},
            {500000,
             {
                 {"alpha", 4.123853571168164},
                 {"beta", 2.07914339855886},
                 {"ct0", 6.052001156627106e-09},
                 {"ct1", 1.3731037996663199e-11},
                 {"ct2", 3.315293060549009e-14},
                 {"k", 9.812805086724211e-08},
             }},
            {800000,
             {
                 {"alpha", 4.123853571168164},
                 {"beta", 2.07914339855886},
                 {"ct0", 6.052001156627106e-09},
                 {"ct1", 1.3731037996663199e-11},
                 {"ct2", 3.315293060549009e-14},
                 {"k", 9.812805086724211e-08},
             }},
            {300000,
             {
                 {"alpha", 2.1466817808584047},
                 {"beta", 2.139630173220604},
                 {"ct0", 0.001541105798542162},
                 {"ct1", 9.25135684683062e-06},
                 {"ct2", 4.791428071078337e-08},
                 {"k", 0.08546283360358237},
             }},
            {400000,
             {
                 {"alpha", 8.27274894957738},
                 {"beta", 2.1617155818557303},
                 {"ct0", 1.3466931115971253e-17},
                 {"ct1", 2.4676217094615232e-20},
                 {"ct2", -3.9377121639760775e-23},
                 {"k", 3.5313702251053173e-22},
             }},
            {50000,
             {
                 {"alpha", 2.3525962926885864},
                 {"beta", 2.5560896822334587},
                 {"ct0", 0.0026688469099427486},
                 {"ct1", 4.061437463160513e-05},
                 {"ct2", 2.1183441402835775e-07},
                 {"k", 0.07070975516481073},
             }},
            {25000,
             {
                 {"alpha", 2.3525962926885864},
                 {"beta", 2.5560896822334587},
                 {"ct0", 0.0026688469099427486},
                 {"ct1", 4.061437463160513e-05},
                 {"ct2", 2.1183441402835775e-07},
                 {"k", 0.07070975516481073},
             }},
            {100000,
             {
                 {"alpha", 1.1716966059510023},
                 {"beta", 2.5226970334652967},
                 {"ct0", 0.003473554573711124},
                 {"ct1", 5.6221032590816286e-05},
                 {"ct2", 3.1740360002324947e-07},
                 {"k", 18956.02913267681},
             }},
            {200000,
             {
                 {"alpha", 3.0647875316227404},
                 {"beta", 2.3232843947192814},
                 {"ct0", 0.0001872613182793334},
                 {"ct1", 1.2193976090778297e-06},
                 {"ct2", 4.743624201506105e-09},
                 {"k", 1.6449252578485438e-05},
             }},
        },
    },
    {
        "N27",
        {
            {1000000,
             {
                 {"alpha", 3.471930172594413},
                 {"beta", 2.1324786840230465},
                 {"ct0", 1.3153259364769057e-07},
                 {"ct1", 5.683033314461967e-13},
                 {"ct2", -3.6567951006427085e-13},
                 {"k", 2.4560291472777438e-05},
             }},
            {500000,
             {
                 {"alpha", 3.4837588474247996},
                 {"beta", 2.1455490665454766},
                 {"ct0", 1.2288678739812055e-07},
                 {"ct1", 6.279227939110073e-11},
                 {"ct2", 1.3557358015564517e-13},
                 {"k", 2.392661765365796e-05},
             }},
            {800000,
             {
                 {"alpha", 3.471930172594413},
                 {"beta", 2.1324786840230465},
                 {"ct0", 1.3153259364769057e-07},
                 {"ct1", 5.683033314461967e-13},
                 {"ct2", -3.6567951006427085e-13},
                 {"k", 2.4560291472777438e-05},
             }},
            {300000,
             {
                 {"alpha", 1.2582003526436814},
                 {"beta", 2.302926598342396},
                 {"ct0", 0.004105630773120781},
                 {"ct1", 3.3646186896806844e-05},
                 {"ct2", 2.128744416819834e-07},
                 {"k", 4253.714257284287},
             }},
            {400000,
             {
                 {"alpha", 3.463391378827239},
                 {"beta", 2.1720528998234827},
                 {"ct0", 5.3316082054133444e-08},
                 {"ct1", -3.41820114808099e-10},
                 {"ct2", -3.3642863526172204e-12},
                 {"k", 7.049993219548801e-05},
             }},
            {50000,
             {
                 {"alpha", 1.529375007546603},
                 {"beta", 2.503764907930501},
                 {"ct0", 0.003020518459475689},
                 {"ct1", 4.8680506435676295e-05},
                 {"ct2", 2.777379127546536e-07},
                 {"k", 592.4756752416615},
             }},
            {25000,
             {
                 {"alpha", 1.529375007546603},
                 {"beta", 2.503764907930501},
                 {"ct0", 0.003020518459475689},
                 {"ct1", 4.8680506435676295e-05},
                 {"ct2", 2.777379127546536e-07},
                 {"k", 592.4756752416615},
             }},
            {100000,
             {
                 {"alpha", 1.9759055367921874},
                 {"beta", 2.5020878353742946},
                 {"ct0", 0.0029204320955145565},
                 {"ct1", 4.6030165427060666e-05},
                 {"ct2", 2.6232471357887897e-07},
                 {"k", 2.957126091831206},
             }},
            {200000,
             {
                 {"alpha", 2.9441963108960847},
                 {"beta", 2.303309838493061},
                 {"ct0", 0.0002449055308139706},
                 {"ct1", 3.2990670154048697e-06},
                 {"ct2", 2.1287286738339177e-08},
                 {"k", 9.04035453836287e-05},
             }},
        },
    },
    {
        "N87",
        {
            {1000000,
             {
                 {"alpha", 4.497128033064765},
                 {"beta", 2.096769106218479},
                 {"ct0", 1.7160388517695716e-09},
                 {"ct1", 8.202553883214351e-12},
                 {"ct2", 8.317657700194287e-14},
                 {"k", 3.815031806374172e-09},
             }},
            {500000,
             {
                 {"alpha", 4.497128033064765},
                 {"beta", 2.096769106218479},
                 {"ct0", 1.7160388517695716e-09},
                 {"ct1", 8.202553883214351e-12},
                 {"ct2", 8.317657700194287e-14},
                 {"k", 3.815031806374172e-09},
             }},
            {800000,
             {
                 {"alpha", 4.497128033064765},
                 {"beta", 2.096769106218479},
                 {"ct0", 1.7160388517695716e-09},
                 {"ct1", 8.202553883214351e-12},
                 {"ct2", 8.317657700194287e-14},
                 {"k", 3.815031806374172e-09},
             }},
            {300000,
             {
                 {"alpha", 1.555208667837359},
                 {"beta", 2.3461121995549883},
                 {"ct0", 0.0030796201235348043},
                 {"ct1", 4.626713131906134e-06},
                 {"ct2", -2.3589772595617744e-08},
                 {"k", 135.28711104911054},
             }},
            {400000,
             {
                 {"alpha", 8.846307238125833},
                 {"beta", 2.1585875625939934},
                 {"ct0", 1.7729462683116735e-18},
                 {"ct1", 1.0110548371210936e-20},
                 {"ct2", 8.998543365427062e-23},
                 {"k", 2.182087218697091e-24},
             }},
            {50000,
             {
                 {"alpha", 2.2583712179625604},
                 {"beta", 2.5443926530330905},
                 {"ct0", 0.002088216273817877},
                 {"ct1", 3.061675126496014e-05},
                 {"ct2", 1.5549853247199972e-07},
                 {"k", 0.4268142896549184},
             }},
            {25000,
             {
                 {"alpha", 2.2583712179625604},
                 {"beta", 2.5443926530330905},
                 {"ct0", 0.002088216273817877},
                 {"ct1", 3.061675126496014e-05},
                 {"ct2", 1.5549853247199972e-07},
                 {"k", 0.4268142896549184},
             }},
            {100000,
             {
                 {"alpha", 1.2688307386196407},
                 {"beta", 2.5541687058102833},
                 {"ct0", 0.003576562713747624},
                 {"ct1", 4.834847003069742e-05},
                 {"ct2", 2.408934879838041e-07},
                 {"k", 9515.310171154813},
             }},
            {200000,
             {
                 {"alpha", 1.3400954633755775},
                 {"beta", 2.385265070709696},
                 {"ct0", 0.0036186538379980805},
                 {"ct1", 6.036145729574879e-05},
                 {"ct2", 4.2011494952919085e-07},
                 {"k", 3244.749891991289},
             }},
        },
    },
    {
        "3C90",
        {
            {1000000,
             {
                 {"alpha", 4.679729716543527},
                 {"beta", 2.108825779644793},
                 {"ct0", 7.482194325920635e-12},
                 {"ct1", -1.910085591357605e-15},
                 {"ct2", 1.232847231127041e-17},
                 {"k", 5.3821928268927866e-08},
             }},
            {500000,
             {
                 {"alpha", 4.679729716543527},
                 {"beta", 2.108825779644793},
                 {"ct0", 7.482194325920635e-12},
                 {"ct1", -1.910085591357605e-15},
                 {"ct2", 1.232847231127041e-17},
                 {"k", 5.3821928268927866e-08},
             }},
            {800000,
             {
                 {"alpha", 4.679729716543527},
                 {"beta", 2.108825779644793},
                 {"ct0", 7.482194325920635e-12},
                 {"ct1", -1.910085591357605e-15},
                 {"ct2", 1.232847231127041e-17},
                 {"k", 5.3821928268927866e-08},
             }},
            {300000,
             {
                 {"alpha", 1.6306179162070384},
                 {"beta", 2.3234439572816745},
                 {"ct0", 0.005426421633737356},
                 {"ct1", 4.511904421784244e-05},
                 {"ct2", 3.2794295460465175e-07},
                 {"k", 21.99141195107354},
             }},
            {400000,
             {
                 {"alpha", 9.292790782848037},
                 {"beta", 2.138804680452066},
                 {"ct0", 2.351520305801877e-13},
                 {"ct1", 9.55159799662906e-16},
                 {"ct2", 8.235603586858213e-18},
                 {"k", 3.2376401933438487e-32},
             }},
            {50000,
             {
                 {"alpha", 1.8569861229888762},
                 {"beta", 2.698697327435223},
                 {"ct0", 0.0029560815064933505},
                 {"ct1", 4.142378323890748e-05},
                 {"ct2", 2.0853044185298125e-07},
                 {"k", 21.500986228654913},
             }},
            {25000,
             {
                 {"alpha", 1.8569861229888762},
                 {"beta", 2.698697327435223},
                 {"ct0", 0.0029560815064933505},
                 {"ct1", 4.142378323890748e-05},
                 {"ct2", 2.0853044185298125e-07},
                 {"k", 21.500986228654913},
             }},
            {100000,
             {
                 {"alpha", 1.4678322191418267},
                 {"beta", 2.67649377134966},
                 {"ct0", 0.0027603633955163753},
                 {"ct1", 4.08025573530478e-05},
                 {"ct2", 2.1347139080662408e-07},
                 {"k", 1191.1729589953318},
             }},
            {200000,
             {
                 {"alpha", 3.480678489960086},
                 {"beta", 2.479282185248303},
                 {"ct0", 1.0403338055430112e-06},
                 {"ct1", -2.476713404645899e-09},
                 {"ct2", -6.670574447108982e-11},
                 {"k", 2.0838234524705946e-05},
             }},
        },
    },
};
std::map<OpenMagnetics::CoreLossesModels, double> maximumAdmittedErrorVolumetricCoreLosses = {
    {OpenMagnetics::CoreLossesModels::STEINMETZ, 1.6}, {OpenMagnetics::CoreLossesModels::IGSE, 1.6},
    {OpenMagnetics::CoreLossesModels::ALBACH, 1.54},   {OpenMagnetics::CoreLossesModels::BARG, 2.1},
    {OpenMagnetics::CoreLossesModels::ROSHEN, 2.48},   {OpenMagnetics::CoreLossesModels::NSE, 1.55},
    {OpenMagnetics::CoreLossesModels::MSE, 1.54},
};
std::map<OpenMagnetics::CoreLossesModels, std::vector<double>> testAverageErrors = {
    {OpenMagnetics::CoreLossesModels::STEINMETZ, {}},
    {OpenMagnetics::CoreLossesModels::IGSE, {}},
    {OpenMagnetics::CoreLossesModels::ALBACH, {}},
};
std::map<OpenMagnetics::CoreLossesModels, double> testMaximumErrors = {
    {OpenMagnetics::CoreLossesModels::STEINMETZ, 0.0},
    {OpenMagnetics::CoreLossesModels::IGSE, 0.0},
    {OpenMagnetics::CoreLossesModels::ALBACH, 0.0},
};


double run_test_core_losses(const OpenMagnetics::CoreLossesModels& modelName,
                            const std::string& shapeName,
                            const std::string& materialName,
                            double frequency,
                            double magneticFluxDensityPeak,
                            double magneticFluxDensityDutyCycle,
                            double magneticFieldStrengthDc,
                            OpenMagnetics::WaveformLabel waveformShape,
                            double temperature,
                            double expectedVolumetricLosses,
                            json steinmetzCoefficients = json({})) {

    auto settings = OpenMagnetics::Settings::GetInstance();
    settings->reset();
    OpenMagnetics::clear_databases();

    double maximumAdmittedErrorVolumetricCoreLossesValue = maximumAdmittedErrorVolumetricCoreLosses[modelName];
    OpenMagnetics::CoreWrapper core = OpenMagneticsTesting::get_quick_core(shapeName, json::array(), 1, materialName);
    auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(modelName);
    if (!steinmetzCoefficients.empty()) {
        OpenMagnetics::SteinmetzCoreLossesMethodRangeDatum steinmetzDatum(steinmetzCoefficients);
        coreLossesModel->set_steinmetz_datum(steinmetzDatum);
    }

    json excitationJson = json();
    excitationJson["frequency"] = frequency;
    excitationJson["magneticFluxDensity"]["processed"]["dutyCycle"] = magneticFluxDensityDutyCycle;
    excitationJson["magneticFluxDensity"]["processed"]["label"] = waveformShape;
    excitationJson["magneticFluxDensity"]["processed"]["offset"] = 0;
    excitationJson["magneticFluxDensity"]["processed"]["peak"] = magneticFluxDensityPeak;
    excitationJson["magneticFluxDensity"]["processed"]["peakToPeak"] = magneticFluxDensityPeak * 2;
    excitationJson["magneticFieldStrength"]["processed"]["offset"] = magneticFieldStrengthDc;
    excitationJson["magneticFieldStrength"]["processed"]["label"] = waveformShape;
    excitationJson["magneticFieldStrength"]["processed"]["peakToPeak"] = 0;

    OpenMagnetics::OperatingPointExcitation excitation(excitationJson);

    auto coreLosses = coreLossesModel->get_core_losses(core, excitation, temperature);
    double calculatedVolumetricCoreLosses = coreLosses.get_volumetric_losses().value();

    double error = abs(expectedVolumetricLosses - calculatedVolumetricCoreLosses) / expectedVolumetricLosses;

    if (error > maximumAdmittedErrorVolumetricCoreLossesValue) {
        // std::cout << "error " << error * 100 << " %" << std::endl;
        // OpenMagneticsTesting::print_json(excitationJson);
    }
    if (testMaximumErrors[modelName] < error) {
        testMaximumErrors[modelName] = error;
    }
    CHECK_CLOSE(calculatedVolumetricCoreLosses, expectedVolumetricLosses,
                expectedVolumetricLosses * maximumAdmittedErrorVolumetricCoreLossesValue);

    return error;
}

void test_core_losses_magnet_verification_3F4(OpenMagnetics::CoreLossesModels modelName,
                                              bool useDynamicCoefficients = false) {
    double meanError = 0;
    std::vector<std::map<std::string, double>> tests = {
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 25000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 35000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 400000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 820000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.05},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 25000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 25000},
            {"magneticFluxDensityPeak", 0.2},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 300000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 1000000},
            {"magneticFluxDensityPeak", 0.02},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 50000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 1000000},
            {"magneticFluxDensityPeak", 0.03},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 150000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 500000},
            {"magneticFluxDensityPeak", 0.05},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 80},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 170000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 500000},
            {"magneticFluxDensityPeak", 0.05},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 50},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 175000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 1000000},
            {"magneticFluxDensityPeak", 0.05},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 600000.0},
        },
    };

    std::string coreShape = "PQ 20/20"; // Provisionally because toroid are not implemented
    std::string coreMaterial = "3F4";
    json steinmetzCoefficients = json({});

    for (auto& test : tests) {
        if (useDynamicCoefficients) {
            steinmetzCoefficients = json(dynamicCoefficients[coreMaterial][test["frequency"]]);
        }
        meanError +=
            run_test_core_losses(modelName, coreShape, coreMaterial, test["frequency"], test["magneticFluxDensityPeak"],
                                 test["magneticFluxDensityDutyCycle"], test["magneticFieldStrengthDc"],
                                 magic_enum::enum_cast<OpenMagnetics::WaveformLabel>(test["waveformShape"]).value(),
                                 test["temperature"], test["expectedVolumetricLosses"], steinmetzCoefficients);
    }
    meanError /= tests.size();
    testAverageErrors[modelName].push_back(meanError);
    if (verboseTests) {
        std::cout << "Mean Error in Core losses for " << coreMaterial << " with Model " << magic_enum::enum_name(modelName)
                  << ": " << meanError * 100 << " %" << std::endl;
        std::cout << "Current average for  " << magic_enum::enum_name(modelName) << ": "
                  << std::reduce(testAverageErrors[modelName].begin(), testAverageErrors[modelName].end()) /
                         testAverageErrors[modelName].size() * 100
                  << " %" << std::endl;
        std::cout << "Current maximum for  " << magic_enum::enum_name(modelName) << ": "
                  << testMaximumErrors[modelName] * 100 << " %" << std::endl;
    }
}

void test_core_losses_magnet_verification_N49(OpenMagnetics::CoreLossesModels modelName,
                                              bool useDynamicCoefficients = false) {
    double meanError = 0;
    std::vector<std::map<std::string, double>> tests = {
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.025},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 2000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.05},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 15000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 100000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 200000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 200000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 200000},
            {"magneticFluxDensityPeak", 0.2},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 1500000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 50000},
            {"magneticFluxDensityPeak", 0.05},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 7000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::TRIANGULAR)},
            {"frequency", 200000},
            {"magneticFluxDensityPeak", 0.15},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 705000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::TRIANGULAR)},
            {"frequency", 50000},
            {"magneticFluxDensityPeak", 0.29},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 1000000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 500000},
            {"magneticFluxDensityPeak", 0.05},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 195650.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 500000},
            {"magneticFluxDensityPeak", 0.05},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 149920.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 500000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 700000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 500000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 80},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 700000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 500000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 550000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 500000},
            {"magneticFluxDensityPeak", 0.0125},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 80},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 2000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 500000},
            {"magneticFluxDensityPeak", 0.025},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 80},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 11000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 500000},
            {"magneticFluxDensityPeak", 0.05},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 80},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 65000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 500000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 80},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 550000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 50000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 44000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 100000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 200000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 220000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 400000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 500000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 800000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 1020000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::TRIANGULAR)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.15},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityPeakOffset", 0.15},
            {"magneticFieldStrengthDc", 20},
            {"expectedVolumetricLosses", 3500000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::TRIANGULAR)},
            {"frequency", 50000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityPeakOffset", 0.15},
            {"magneticFieldStrengthDc", 20},
            {"expectedVolumetricLosses", 700000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::TRIANGULAR)},
            {"frequency", 400000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityPeakOffset", 0.2},
            {"magneticFieldStrengthDc", 25},
            {"expectedVolumetricLosses", 9500005.0},
        },
    };

    std::string coreShape = "PQ 20/20"; // Provisionally because toroid are not implemented
    std::string coreMaterial = "N49";
    json steinmetzCoefficients = json({});

    for (auto& test : tests) {
        if (useDynamicCoefficients) {
            steinmetzCoefficients = json(dynamicCoefficients[coreMaterial][test["frequency"]]);
        }
        meanError +=
            run_test_core_losses(modelName, coreShape, coreMaterial, test["frequency"], test["magneticFluxDensityPeak"],
                                 test["magneticFluxDensityDutyCycle"], test["magneticFieldStrengthDc"],
                                 magic_enum::enum_cast<OpenMagnetics::WaveformLabel>(test["waveformShape"]).value(),
                                 test["temperature"], test["expectedVolumetricLosses"], steinmetzCoefficients);
    }
    meanError /= tests.size();
    testAverageErrors[modelName].push_back(meanError);
    if (verboseTests) {
        std::cout << "Mean Error in Core losses for " << coreMaterial << " with Model " << magic_enum::enum_name(modelName)
                  << ": " << meanError * 100 << " %" << std::endl;
        std::cout << "Current average for  " << magic_enum::enum_name(modelName) << ": "
                  << std::reduce(testAverageErrors[modelName].begin(), testAverageErrors[modelName].end()) /
                         testAverageErrors[modelName].size() * 100
                  << " %" << std::endl;
        std::cout << "Current maximum for  " << magic_enum::enum_name(modelName) << ": "
                  << testMaximumErrors[modelName] * 100 << " %" << std::endl;
    }
}

void test_core_losses_magnet_verification_3C94(OpenMagnetics::CoreLossesModels modelName,
                                               bool useDynamicCoefficients = false) {
    double meanError = 0;
    std::vector<std::map<std::string, double>> tests = {
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.2},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 300000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 200000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 160000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 300000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 1050000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 200000},
            {"magneticFluxDensityPeak", 0.2},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 1050000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 25000},
            {"magneticFluxDensityPeak", 0.2},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 60000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 25000},
            {"magneticFluxDensityPeak", 0.3},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 190000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.08},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 30000.0},
        },
    };

    std::string coreShape = "PQ 20/20"; // Provisionally because toroid are not implemented
    std::string coreMaterial = "3C94";
    json steinmetzCoefficients = json({});

    for (auto& test : tests) {
        if (useDynamicCoefficients) {
            steinmetzCoefficients = json(dynamicCoefficients[coreMaterial][test["frequency"]]);
        }
        meanError +=
            run_test_core_losses(modelName, coreShape, coreMaterial, test["frequency"], test["magneticFluxDensityPeak"],
                                 test["magneticFluxDensityDutyCycle"], test["magneticFieldStrengthDc"],
                                 magic_enum::enum_cast<OpenMagnetics::WaveformLabel>(test["waveformShape"]).value(),
                                 test["temperature"], test["expectedVolumetricLosses"], steinmetzCoefficients);
    }
    meanError /= tests.size();
    testAverageErrors[modelName].push_back(meanError);
    if (verboseTests) {
        std::cout << "Mean Error in Core losses for " << coreMaterial << " with Model " << magic_enum::enum_name(modelName)
                  << ": " << meanError * 100 << " %" << std::endl;
        std::cout << "Current average for  " << magic_enum::enum_name(modelName) << ": "
                  << std::reduce(testAverageErrors[modelName].begin(), testAverageErrors[modelName].end()) /
                         testAverageErrors[modelName].size() * 100
                  << " %" << std::endl;
        std::cout << "Current maximum for  " << magic_enum::enum_name(modelName) << ": "
                  << testMaximumErrors[modelName] * 100 << " %" << std::endl;
    }
}

void test_core_losses_magnet_verification_N27(OpenMagnetics::CoreLossesModels modelName,
                                              bool useDynamicCoefficients = false) {
    double meanError = 0;
    std::vector<std::map<std::string, double>> tests = {
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 25000},
            {"magneticFluxDensityPeak", 0.2},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 155000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.2},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 920000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 25000},
            {"magneticFluxDensityPeak", 0.05},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 60},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 7000.0},
        },
    };

    std::string coreShape = "PQ 20/20"; // Provisionally because toroid are not implemented
    std::string coreMaterial = "N27";
    json steinmetzCoefficients = json({});

    for (auto& test : tests) {
        if (useDynamicCoefficients) {
            steinmetzCoefficients = json(dynamicCoefficients[coreMaterial][test["frequency"]]);
        }
        meanError +=
            run_test_core_losses(modelName, coreShape, coreMaterial, test["frequency"], test["magneticFluxDensityPeak"],
                                 test["magneticFluxDensityDutyCycle"], test["magneticFieldStrengthDc"],
                                 magic_enum::enum_cast<OpenMagnetics::WaveformLabel>(test["waveformShape"]).value(),
                                 test["temperature"], test["expectedVolumetricLosses"], steinmetzCoefficients);
    }
    meanError /= tests.size();
    testAverageErrors[modelName].push_back(meanError);
    if (verboseTests) {
        std::cout << "Mean Error in Core losses for " << coreMaterial << " with Model " << magic_enum::enum_name(modelName)
                  << ": " << meanError * 100 << " %" << std::endl;
        std::cout << "Current average for  " << magic_enum::enum_name(modelName) << ": "
                  << std::reduce(testAverageErrors[modelName].begin(), testAverageErrors[modelName].end()) /
                         testAverageErrors[modelName].size() * 100
                  << " %" << std::endl;
        std::cout << "Current maximum for  " << magic_enum::enum_name(modelName) << ": "
                  << testMaximumErrors[modelName] * 100 << " %" << std::endl;
    }
}

void test_core_losses_magnet_verification_N87(OpenMagnetics::CoreLossesModels modelName,
                                              bool useDynamicCoefficients = false) {
    double meanError = 0;
    std::vector<std::map<std::string, double>> tests = {
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::TRIANGULAR)},
            {"frequency", 50000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 62000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::TRIANGULAR)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.24},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 1000000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::TRIANGULAR)},
            {"frequency", 400000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 900000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::TRIANGULAR)},
            {"frequency", 400000},
            {"magneticFluxDensityPeak", 0.05},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 180000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.025},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 1500.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.05},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 40},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 20000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.05},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 90},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 8000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 80},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 60000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.2},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 30},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 600000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::TRIANGULAR)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.2},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityPeakOffset", 0.2},
            {"magneticFieldStrengthDc", 30},
            {"expectedVolumetricLosses", 9000000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::TRIANGULAR)},
            {"frequency", 50000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityPeakOffset", 0.2},
            {"magneticFieldStrengthDc", 30},
            {"expectedVolumetricLosses", 1000000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.2},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 608400.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.2},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 622660.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::TRIANGULAR)},
            {"frequency", 400000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityPeakOffset", 0.3},
            {"magneticFieldStrengthDc", 35},
            {"expectedVolumetricLosses", 11700005.0},
        },
    };

    std::string coreShape = "PQ 20/20"; // Provisionally because toroid are not implemented
    std::string coreMaterial = "N87";
    json steinmetzCoefficients = json({});

    for (auto& test : tests) {
        if (useDynamicCoefficients) {
            steinmetzCoefficients = json(dynamicCoefficients[coreMaterial][test["frequency"]]);
        }
        meanError +=
            run_test_core_losses(modelName, coreShape, coreMaterial, test["frequency"], test["magneticFluxDensityPeak"],
                                 test["magneticFluxDensityDutyCycle"], test["magneticFieldStrengthDc"],
                                 magic_enum::enum_cast<OpenMagnetics::WaveformLabel>(test["waveformShape"]).value(),
                                 test["temperature"], test["expectedVolumetricLosses"], steinmetzCoefficients);
    }
    meanError /= tests.size();
    testAverageErrors[modelName].push_back(meanError);
    if (verboseTests) {
        std::cout << "Mean Error in Core losses for " << coreMaterial << " with Model " << magic_enum::enum_name(modelName)
                  << ": " << meanError * 100 << " %" << std::endl;
        std::cout << "Current average for  " << magic_enum::enum_name(modelName) << ": "
                  << std::reduce(testAverageErrors[modelName].begin(), testAverageErrors[modelName].end()) /
                         testAverageErrors[modelName].size() * 100
                  << " %" << std::endl;
        std::cout << "Current maximum for  " << magic_enum::enum_name(modelName) << ": "
                  << testMaximumErrors[modelName] * 100 << " %" << std::endl;
    }
}

void test_core_losses_magnet_verification_3C90(OpenMagnetics::CoreLossesModels modelName,
                                               bool useDynamicCoefficients = false) {
    double meanError = 0;
    std::vector<std::map<std::string, double>> tests = {
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::TRIANGULAR)},
            {"frequency", 50000},
            {"magneticFluxDensityPeak", 0.05},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 10000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::TRIANGULAR)},
            {"frequency", 400000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 25},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 895000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 99530.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 25000},
            {"magneticFluxDensityPeak", 0.2},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 121000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 111670.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 200000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 250000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 50000},
            {"magneticFluxDensityPeak", 0.2},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 200000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.2},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 480000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.2},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 480000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 25000},
            {"magneticFluxDensityPeak", 0.3},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 200000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 100000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 100},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 70000.0},
        },
        {
            {"waveformShape", magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL)},
            {"frequency", 200000},
            {"magneticFluxDensityPeak", 0.1},
            {"magneticFluxDensityDutyCycle", 0.5},
            {"temperature", 60},
            {"magneticFluxDensityOffset", 0},
            {"magneticFieldStrengthDc", 0},
            {"expectedVolumetricLosses", 300000.0},
        },
    };

    std::string coreShape = "PQ 20/20"; // Provisionally because toroid are not implemented
    std::string coreMaterial = "3C90";
    json steinmetzCoefficients = json({});

    for (auto& test : tests) {
        if (useDynamicCoefficients) {
            steinmetzCoefficients = json(dynamicCoefficients[coreMaterial][test["frequency"]]);
        }
        meanError +=
            run_test_core_losses(modelName, coreShape, coreMaterial, test["frequency"], test["magneticFluxDensityPeak"],
                                 test["magneticFluxDensityDutyCycle"], test["magneticFieldStrengthDc"],
                                 magic_enum::enum_cast<OpenMagnetics::WaveformLabel>(test["waveformShape"]).value(),
                                 test["temperature"], test["expectedVolumetricLosses"], steinmetzCoefficients);
    }
    meanError /= tests.size();
    testAverageErrors[modelName].push_back(meanError);
    if (verboseTests) {
        std::cout << "Mean Error in Core losses for " << coreMaterial << " with Model " << magic_enum::enum_name(modelName)
                  << ": " << meanError * 100 << " %" << std::endl;
        std::cout << "Current average for  " << magic_enum::enum_name(modelName) << ": "
                  << std::reduce(testAverageErrors[modelName].begin(), testAverageErrors[modelName].end()) /
                         testAverageErrors[modelName].size() * 100
                  << " %" << std::endl;
        std::cout << "Current maximum for  " << magic_enum::enum_name(modelName) << ": "
                  << testMaximumErrors[modelName] * 100 << " %" << std::endl;
    }
}

std::vector<std::map<std::string, double>> load_sample_data_from_material(std::string material) {
    std::string file_path = __FILE__;
    auto sample_path = file_path.substr(0, file_path.rfind("/")).append("/testData/" + material + "_sample.csv");
    // auto sample_path = file_path.substr(0, file_path.rfind("/")).append("/testData/" + material + "_database.csv");

    std::ifstream in(sample_path);
    std::vector<std::vector<double>> fields;
    bool headers = false;
    size_t number_read_rows = 0;

    if (in) {
        std::string line;

        while (getline(in, line)) {
            std::stringstream sep(line);
            std::string field;

            if (!headers) {
                headers = true;
                continue;
            }

            std::vector<double> row_data;

            while (getline(sep, field, ',')) {
                row_data.push_back(stod(field));
            }

            fields.push_back(row_data);
            number_read_rows++;
            // if (number_read_rows > 500)
            // break;
        }
    }

    std::vector<std::map<std::string, double>> tests;
    for (auto row : fields) {
        std::map<std::string, double> elem;
        double dutyCycle;
        if (row[3] == -1) {
            elem["waveformShape"] = magic_enum::enum_integer(OpenMagnetics::WaveformLabel::SINUSOIDAL);
            dutyCycle = 0.5;
        }
        else if (row[3] + row[4] == 1) {
            elem["waveformShape"] = magic_enum::enum_integer(OpenMagnetics::WaveformLabel::TRIANGULAR);
            dutyCycle = row[3];
        }
        else {
            continue;
        }

        elem["frequency"] = row[0];
        elem["magneticFluxDensityPeak"] = row[1];
        elem["magneticFluxDensityDutyCycle"] = dutyCycle;
        elem["temperature"] = row[5];
        elem["magneticFieldStrengthDc"] = row[2];
        elem["expectedVolumetricLosses"] = row[6];

        tests.push_back(elem);
    }
    return tests;
}

void export_test_result_for_material(std::vector<std::map<std::string, double>> tests,
                                     std::string material,
                                     OpenMagnetics::CoreLossesModels modelName) {
    std::string file_path = __FILE__;
    auto results_path =
        file_path.substr(0, file_path.rfind("/"))
            .append("/testData/" + std::string(magic_enum::enum_name(modelName)) + "_" + material + "_result.csv");

    std::vector<std::vector<double>> fields;

    std::vector<std::string> headers;
    for (auto const& element : tests[0]) {
        headers.push_back(element.first);
    }
    std::ofstream myfile;
    myfile.open(results_path);
    for (auto& header : headers) {
        myfile << header;
        myfile << ",";
    }
    myfile << "\n";

    for (auto const& test : tests) {
        for (auto const& element : test) {
            myfile << element.second;
            myfile << ",";
        }
        myfile << "\n";
    }
    myfile.close();
}

void test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels modelName,
                                  std::string coreMaterial,
                                  bool useDynamicCoefficients = false) {
    double meanError = 0;
    std::string coreShape = "PQ 20/20"; // Provisionally because toroid are not implemented
    auto tests = load_sample_data_from_material(coreMaterial);
    std::vector<std::map<std::string, double>> testResult;

    json steinmetzCoefficients = json({});
    for (auto& test : tests) {
        // OpenMagneticsTesting::print(test["expectedVolumetricLosses"], steinmetzCoefficients);

        if (useDynamicCoefficients) {
            steinmetzCoefficients = json(dynamicCoefficients[coreMaterial][test["frequency"]]);
        }
        double testError =
            run_test_core_losses(modelName, coreShape, coreMaterial, test["frequency"], test["magneticFluxDensityPeak"],
                                 test["magneticFluxDensityDutyCycle"], test["magneticFieldStrengthDc"],
                                 magic_enum::enum_cast<OpenMagnetics::WaveformLabel>(test["waveformShape"]).value(),
                                 test["temperature"], test["expectedVolumetricLosses"], steinmetzCoefficients);
        test["error"] = testError;
        testResult.push_back(test);
        meanError += testError;
    }
    meanError /= tests.size();
    testAverageErrors[modelName].push_back(meanError);
    if (verboseTests) {
        std::cout << "Mean Error in Core losses for " << coreMaterial << " with Model " << magic_enum::enum_name(modelName)
                  << ": " << meanError * 100 << " %" << std::endl;
        std::cout << "Current average for  " << magic_enum::enum_name(modelName) << ": "
                  << std::reduce(testAverageErrors[modelName].begin(), testAverageErrors[modelName].end()) /
                         testAverageErrors[modelName].size() * 100
                  << " %" << std::endl;
        std::cout << "Current maximum for  " << magic_enum::enum_name(modelName) << ": "
                  << testMaximumErrors[modelName] * 100 << " %" << std::endl;
    }

    export_test_result_for_material(testResult, coreMaterial, modelName);
}

SUITE(SteinmetzModel) {
    TEST(Test_PQ_20_20_3F4) {
        test_core_losses_magnet_verification_3F4(OpenMagnetics::CoreLossesModels::STEINMETZ);
    }
    TEST(Test_PQ_20_20_N49) {
        test_core_losses_magnet_verification_N49(OpenMagnetics::CoreLossesModels::STEINMETZ);
    }
    TEST(Test_PQ_20_20_3C94) {
        test_core_losses_magnet_verification_3C94(OpenMagnetics::CoreLossesModels::STEINMETZ);
    }
    TEST(Test_PQ_20_20_N27) {
        test_core_losses_magnet_verification_N27(OpenMagnetics::CoreLossesModels::STEINMETZ);
    }
    TEST(Test_PQ_20_20_N87) {
        test_core_losses_magnet_verification_N87(OpenMagnetics::CoreLossesModels::STEINMETZ);
    }
    TEST(Test_PQ_20_20_3C90) {
        test_core_losses_magnet_verification_3C90(OpenMagnetics::CoreLossesModels::STEINMETZ);
    }
    TEST(Test_Magnet_3C90) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::STEINMETZ, "3C90");
    }
    TEST(Test_Magnet_3C94) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::STEINMETZ, "3C94");
    }
    TEST(Test_Magnet_3F4) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::STEINMETZ, "3F4");
    }
    TEST(Test_Magnet_N27) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::STEINMETZ, "N27");
    }
    TEST(Test_Magnet_N30) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::STEINMETZ, "N30");
    }
    TEST(Test_Magnet_N49) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::STEINMETZ, "N49");
    }
    TEST(Test_Magnet_N87) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::STEINMETZ, "N87");
    }

    // TEST(dynamic_coefficients)
    // {
    //     test_core_losses_magnet_verification_3F4(OpenMagnetics::CoreLossesModels::STEINMETZ, true);
    //     test_core_losses_magnet_verification_N49(OpenMagnetics::CoreLossesModels::STEINMETZ, true);
    //     test_core_losses_magnet_verification_3C94(OpenMagnetics::CoreLossesModels::STEINMETZ, true);
    //     test_core_losses_magnet_verification_N27(OpenMagnetics::CoreLossesModels::STEINMETZ, true);
    //     test_core_losses_magnet_verification_N87(OpenMagnetics::CoreLossesModels::STEINMETZ, true);
    //     test_core_losses_magnet_verification_3C90(OpenMagnetics::CoreLossesModels::STEINMETZ, true);
    // }
}

SUITE(IGSEModel) {
    TEST(Test_Ki_3C95) {
        std::string shapeName = "PQ 20/20"; // Provisionally because toroid are not implemented
        std::string materialName = "3C95";

        OpenMagnetics::CoreWrapper core = OpenMagneticsTesting::get_quick_core(shapeName, json::array(), 1, materialName);
        auto steinmetzDatum = OpenMagnetics::CoreLossesModel::get_steinmetz_coefficients(
            core.get_functional_description().get_material(), 100000);
        auto coreLossesIGSEModel = OpenMagnetics::CoreLossesIGSEModel();

        auto ki = coreLossesIGSEModel.get_ki(steinmetzDatum);
        double expectedKi = 8.17;

        CHECK_CLOSE(ki, expectedKi, expectedKi * 0.1);
    }
    TEST(Test_PQ_20_20_3F4) {
        test_core_losses_magnet_verification_3F4(OpenMagnetics::CoreLossesModels::IGSE);
    }
    TEST(Test_PQ_20_20_N49) {
        test_core_losses_magnet_verification_N49(OpenMagnetics::CoreLossesModels::IGSE);
    }
    TEST(Test_PQ_20_20_3C94) {
        test_core_losses_magnet_verification_3C94(OpenMagnetics::CoreLossesModels::IGSE);
    }
    TEST(Test_PQ_20_20_N27) {
        test_core_losses_magnet_verification_N27(OpenMagnetics::CoreLossesModels::IGSE);
    }
    TEST(Test_PQ_20_20_N87) {
        test_core_losses_magnet_verification_N87(OpenMagnetics::CoreLossesModels::IGSE);
    }
    TEST(Test_PQ_20_20_3C90) {
        test_core_losses_magnet_verification_3C90(OpenMagnetics::CoreLossesModels::IGSE);
    }
    TEST(Test_Magnet_3C90) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::IGSE, "3C90");
    }
    TEST(Test_Magnet_3C94) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::IGSE, "3C94");
    }
    TEST(Test_Magnet_3F4) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::IGSE, "3F4");
    }
    TEST(Test_Magnet_N27) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::IGSE, "N27");
    }
    TEST(Test_Magnet_N30) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::IGSE, "N30");
    }
    TEST(Test_Magnet_N49) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::IGSE, "N49");
    }
    TEST(Test_Magnet_N87) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::IGSE, "N87");
    }
    // TEST(dynamic_coefficients)
    // {
    //     test_core_losses_magnet_verification_3F4(OpenMagnetics::CoreLossesModels::IGSE, true);
    //     test_core_losses_magnet_verification_N49(OpenMagnetics::CoreLossesModels::IGSE, true);
    //     test_core_losses_magnet_verification_3C94(OpenMagnetics::CoreLossesModels::IGSE, true);
    //     test_core_losses_magnet_verification_N27(OpenMagnetics::CoreLossesModels::IGSE, true);
    //     test_core_losses_magnet_verification_N87(OpenMagnetics::CoreLossesModels::IGSE, true);
    //     test_core_losses_magnet_verification_3C90(OpenMagnetics::CoreLossesModels::IGSE, true);
    // }
}

SUITE(AlbachModel) {
    TEST(Test_PQ_20_20_3F4) {
        test_core_losses_magnet_verification_3F4(OpenMagnetics::CoreLossesModels::ALBACH);
    }
    TEST(Test_PQ_20_20_N49) {
        test_core_losses_magnet_verification_N49(OpenMagnetics::CoreLossesModels::ALBACH);
    }
    TEST(Test_PQ_20_20_3C94) {
        test_core_losses_magnet_verification_3C94(OpenMagnetics::CoreLossesModels::ALBACH);
    }
    TEST(Test_PQ_20_20_N27) {
        test_core_losses_magnet_verification_N27(OpenMagnetics::CoreLossesModels::ALBACH);
    }
    TEST(Test_PQ_20_20_N87) {
        test_core_losses_magnet_verification_N87(OpenMagnetics::CoreLossesModels::ALBACH);
    }
    TEST(Test_PQ_20_20_3C90) {
        test_core_losses_magnet_verification_3C90(OpenMagnetics::CoreLossesModels::ALBACH);
    }
    TEST(Test_Magnet_3C90) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::ALBACH, "3C90");
    }
    TEST(Test_Magnet_3C94) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::ALBACH, "3C94");
    }
    TEST(Test_Magnet_3F4) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::ALBACH, "3F4");
    }
    TEST(Test_Magnet_N27) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::ALBACH, "N27");
    }
    TEST(Test_Magnet_N30) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::ALBACH, "N30");
    }
    TEST(Test_Magnet_N49) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::ALBACH, "N49");
    }
    TEST(Test_Magnet_N87) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::ALBACH, "N87");
    }
    // TEST(dynamic_coefficients)
    // {
    //     test_core_losses_magnet_verification_3F4(OpenMagnetics::CoreLossesModels::ALBACH, true);
    //     test_core_losses_magnet_verification_N49(OpenMagnetics::CoreLossesModels::ALBACH, true);
    //     test_core_losses_magnet_verification_3C94(OpenMagnetics::CoreLossesModels::ALBACH, true);
    //     test_core_losses_magnet_verification_N27(OpenMagnetics::CoreLossesModels::ALBACH, true);
    //     test_core_losses_magnet_verification_N87(OpenMagnetics::CoreLossesModels::ALBACH, true);
    //     test_core_losses_magnet_verification_3C90(OpenMagnetics::CoreLossesModels::ALBACH, true);
    // }
}

SUITE(MSEModel) {
    TEST(Test_PQ_20_20_3F4) {
        test_core_losses_magnet_verification_3F4(OpenMagnetics::CoreLossesModels::MSE);
    }
    TEST(Test_PQ_20_20_N49) {
        test_core_losses_magnet_verification_N49(OpenMagnetics::CoreLossesModels::MSE);
    }
    TEST(Test_PQ_20_20_3C94) {
        test_core_losses_magnet_verification_3C94(OpenMagnetics::CoreLossesModels::MSE);
    }
    TEST(Test_PQ_20_20_N27) {
        test_core_losses_magnet_verification_N27(OpenMagnetics::CoreLossesModels::MSE);
    }
    TEST(Test_PQ_20_20_N87) {
        test_core_losses_magnet_verification_N87(OpenMagnetics::CoreLossesModels::MSE);
    }
    TEST(Test_PQ_20_20_3C90) {
        test_core_losses_magnet_verification_3C90(OpenMagnetics::CoreLossesModels::MSE);
    }
    TEST(Test_Magnet_3C90) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::MSE, "3C90");
    }
    TEST(Test_Magnet_3C94) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::MSE, "3C94");
    }
    TEST(Test_Magnet_3F4) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::MSE, "3F4");
    }
    TEST(Test_Magnet_N27) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::MSE, "N27");
    }
    TEST(Test_Magnet_N30) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::MSE, "N30");
    }
    TEST(Test_Magnet_N49) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::MSE, "N49");
    }
    TEST(Test_Magnet_N87) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::MSE, "N87");
    }
    // TEST(dynamic_coefficients)
    // {
    //     test_core_losses_magnet_verification_3F4(OpenMagnetics::CoreLossesModels::MSE, true);
    //     test_core_losses_magnet_verification_N49(OpenMagnetics::CoreLossesModels::MSE, true);
    //     test_core_losses_magnet_verification_3C94(OpenMagnetics::CoreLossesModels::MSE, true);
    //     test_core_losses_magnet_verification_N27(OpenMagnetics::CoreLossesModels::MSE, true);
    //     test_core_losses_magnet_verification_N87(OpenMagnetics::CoreLossesModels::MSE, true);
    //     test_core_losses_magnet_verification_3C90(OpenMagnetics::CoreLossesModels::MSE, true);
    // }
}

SUITE(NSEModel) {
    TEST(Test_PQ_20_20_3F4) {
        test_core_losses_magnet_verification_3F4(OpenMagnetics::CoreLossesModels::NSE);
    }
    TEST(Test_PQ_20_20_N49) {
        test_core_losses_magnet_verification_N49(OpenMagnetics::CoreLossesModels::NSE);
    }
    TEST(Test_PQ_20_20_3C94) {
        test_core_losses_magnet_verification_3C94(OpenMagnetics::CoreLossesModels::NSE);
    }
    TEST(Test_PQ_20_20_N27) {
        test_core_losses_magnet_verification_N27(OpenMagnetics::CoreLossesModels::NSE);
    }
    TEST(Test_PQ_20_20_N87) {
        test_core_losses_magnet_verification_N87(OpenMagnetics::CoreLossesModels::NSE);
    }
    TEST(Test_PQ_20_20_3C90) {
        test_core_losses_magnet_verification_3C90(OpenMagnetics::CoreLossesModels::NSE);
    }
    TEST(Test_Magnet_3C90) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::NSE, "3C90");
    }
    TEST(Test_Magnet_3C94) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::NSE, "3C94");
    }
    TEST(Test_Magnet_3F4) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::NSE, "3F4");
    }
    TEST(Test_Magnet_N27) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::NSE, "N27");
    }
    TEST(Test_Magnet_N30) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::NSE, "N30");
    }
    TEST(Test_Magnet_N49) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::NSE, "N49");
    }
    TEST(Test_Magnet_N87) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::NSE, "N87");
    }
    // TEST(dynamic_coefficients)
    // {
    //     test_core_losses_magnet_verification_3F4(OpenMagnetics::CoreLossesModels::NSE, true);
    //     test_core_losses_magnet_verification_N49(OpenMagnetics::CoreLossesModels::NSE, true);
    //     test_core_losses_magnet_verification_3C94(OpenMagnetics::CoreLossesModels::NSE, true);
    //     test_core_losses_magnet_verification_N27(OpenMagnetics::CoreLossesModels::NSE, true);
    //     test_core_losses_magnet_verification_N87(OpenMagnetics::CoreLossesModels::NSE, true);
    //     test_core_losses_magnet_verification_3C90(OpenMagnetics::CoreLossesModels::NSE, true);
    // }
}

SUITE(BargModel) {
    TEST(Test_PQ_20_20_3F4) {
        test_core_losses_magnet_verification_3F4(OpenMagnetics::CoreLossesModels::BARG);
    }
    TEST(Test_PQ_20_20_N49) {
        test_core_losses_magnet_verification_N49(OpenMagnetics::CoreLossesModels::BARG);
    }
    TEST(Test_PQ_20_20_3C94) {
        test_core_losses_magnet_verification_3C94(OpenMagnetics::CoreLossesModels::BARG);
    }
    TEST(Test_PQ_20_20_N27) {
        test_core_losses_magnet_verification_N27(OpenMagnetics::CoreLossesModels::BARG);
    }
    TEST(Test_PQ_20_20_N87) {
        test_core_losses_magnet_verification_N87(OpenMagnetics::CoreLossesModels::BARG);
    }
    TEST(Test_PQ_20_20_3C90) {
        test_core_losses_magnet_verification_3C90(OpenMagnetics::CoreLossesModels::BARG);
    }
    TEST(Test_Magnet_3C90) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::BARG, "3C90");
    }
    TEST(Test_Magnet_3C94) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::BARG, "3C94");
    }
    TEST(Test_Magnet_3F4) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::BARG, "3F4");
    }
    TEST(Test_Magnet_N27) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::BARG, "N27");
    }
    TEST(Test_Magnet_N30) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::BARG, "N30");
    }
    TEST(Test_Magnet_N49) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::BARG, "N49");
    }
    TEST(Test_Magnet_N87) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::BARG, "N87");
    }
    // TEST(dynamic_coefficients)
    // {
    //     test_core_losses_magnet_verification_3F4(OpenMagnetics::CoreLossesModels::BARG, true);
    //     test_core_losses_magnet_verification_N49(OpenMagnetics::CoreLossesModels::BARG, true);
    //     test_core_losses_magnet_verification_3C94(OpenMagnetics::CoreLossesModels::BARG, true);
    //     test_core_losses_magnet_verification_N27(OpenMagnetics::CoreLossesModels::BARG, true);
    //     test_core_losses_magnet_verification_N87(OpenMagnetics::CoreLossesModels::BARG, true);
    //     test_core_losses_magnet_verification_3C90(OpenMagnetics::CoreLossesModels::BARG, true);
    // }
}

SUITE(RoshenModel) {
    TEST(Test_PQ_20_20_N49) {
        test_core_losses_magnet_verification_N49(OpenMagnetics::CoreLossesModels::ROSHEN);
    }
    TEST(Test_PQ_20_20_3C94) {
        test_core_losses_magnet_verification_3C94(OpenMagnetics::CoreLossesModels::ROSHEN);
    }
    TEST(Test_PQ_20_20_N27) {
        test_core_losses_magnet_verification_N27(OpenMagnetics::CoreLossesModels::ROSHEN);
    }
    TEST(Test_PQ_20_20_N87) {
        test_core_losses_magnet_verification_N87(OpenMagnetics::CoreLossesModels::ROSHEN);
    }
    TEST(Test_PQ_20_20_3C90) {
        test_core_losses_magnet_verification_3C90(OpenMagnetics::CoreLossesModels::ROSHEN);
    }
    TEST(Test_Magnet_3C90) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::ROSHEN, "3C90");
    }
    TEST(Test_Magnet_3C94) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::ROSHEN, "3C94");
    }
    TEST(Test_Magnet_N27) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::ROSHEN, "N27");
    }
    TEST(Test_Magnet_N49) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::ROSHEN, "N49");
    }
    TEST(Test_Magnet_N87) {
        test_core_losses_magnet_data(OpenMagnetics::CoreLossesModels::ROSHEN, "N87");
    }
    // TEST(dynamic_coefficients)
    // {
    //     test_core_losses_magnet_verification_N49(OpenMagnetics::CoreLossesModels::ROSHEN, true);
    //     test_core_losses_magnet_verification_3C94(OpenMagnetics::CoreLossesModels::ROSHEN, true);
    //     test_core_losses_magnet_verification_N27(OpenMagnetics::CoreLossesModels::ROSHEN, true);
    //     test_core_losses_magnet_verification_N87(OpenMagnetics::CoreLossesModels::ROSHEN, true);
    //     test_core_losses_magnet_verification_3C90(OpenMagnetics::CoreLossesModels::ROSHEN, true);
    // }
}

SUITE(CoreLossesAssorted) {
    double maxError = 0.05;
    TEST(Voltage_And_Current) {
        auto models = json::parse("{\"coreLosses\": \"IGSE\", \"gapReluctance\": \"BALAKRISHNAN\"}");
        auto core = OpenMagnetics::CoreWrapper(json::parse(
            "{\"functionalDescription\": {\"gapping\": [{\"area\": 9.8e-05, \"coordinates\": [0.0, 0.0001, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011301, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 0.0002, \"sectionDimensions\": [0.0092, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"subtractive\"}, {\"area\": 4.7e-05, \"coordinates\": [0.0138, 0.0, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011498, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.004401, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"residual\"}, {\"area\": 4.7e-05, \"coordinates\": [-0.0138, 0.0, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011498, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.004401, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"residual\"}], \"material\": \"3C97\", \"numberStacks\": 1, \"shape\": {\"aliases\": [], \"dimensions\": "
            "{\"A\": 0.032, \"B\": 0.0161, \"C\": 0.01065, \"D\": 0.0115, \"E\": 0.0232, \"F\": 0.0092, \"G\": 0.0, "
            "\"H\": 0.0}, \"family\": \"e\", \"familySubtype\": null, \"magneticCircuit\": \"open\", \"name\": \"E "
            "32/16/11\", \"type\": \"standard\"}, \"type\": \"two-piece set\"}, \"manufacturerInfo\": null, \"name\": "
            "\"My Core\"}"));
        auto winding = OpenMagnetics::CoilWrapper(
            json::parse("{\"bobbin\": \"Dummy\", \"functionalDescription\": [{\"isolationSide\": \"primary\", \"name\": "
                        "\"Primary\", \"numberParallels\": 1, \"numberTurns\": 33, \"wire\": \"Dummy\"}], "
                        "\"layersDescription\": null, \"sectionsDescription\": null, \"turnsDescription\": null}"));
        auto operatingPoint = OpenMagnetics::OperatingPoint(json::parse(
            "{\"conditions\": {\"ambientRelativeHumidity\": null, \"ambientTemperature\": 37.0, \"cooling\": null, "
            "\"name\": null}, \"excitationsPerWinding\": [{\"current\": {\"harmonics\": null, \"processed\": null, "
            "\"waveform\": {\"ancillaryLabel\": null, \"data\": [-5.0, 5.0, -5.0], \"numberPeriods\": null, \"time\": "
            "[0.0, 2.4999999999999998e-06, 1e-05]}}, \"frequency\": 100000.0, \"magneticFieldStrength\": null, "
            "\"magneticFluxDensity\": null, \"magnetizingCurrent\": null, \"name\": \"My Operating Point\", "
            "\"voltage\": {\"harmonics\": null, \"processed\": null, \"waveform\": {\"ancillaryLabel\": null, "
            "\"data\": [7.5, 7.5, -2.4999999999999996, -2.4999999999999996, 7.5], \"numberPeriods\": null, \"time\": "
            "[0.0, 2.4999999999999998e-06, 2.4999999999999998e-06, 1e-05, 1e-05]}}}], \"name\": null}"));

        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::string{models["gapReluctance"]});

        OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];

        auto magneticFluxDensity =
            magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint).second;

        excitation.set_magnetic_flux_density(magneticFluxDensity);
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();

        auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
        auto coreLosses = coreLossesModel->get_core_losses(core, excitation, temperature);
    }
    TEST(Only_Current_0) {
        auto models = json::parse("{\"coreLosses\": \"IGSE\", \"gapReluctance\": \"BALAKRISHNAN\"}");
        auto core = OpenMagnetics::CoreWrapper(json::parse(
            "{\"functionalDescription\": {\"gapping\": [{\"area\": 9.8e-05, \"coordinates\": [0.0, 0.0001, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011301, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 0.0002, \"sectionDimensions\": [0.0092, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"subtractive\"}, {\"area\": 4.7e-05, \"coordinates\": [0.0138, 0.0, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011498, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.004401, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"residual\"}, {\"area\": 4.7e-05, \"coordinates\": [-0.0138, 0.0, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011498, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.004401, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"residual\"}], \"material\": \"3C97\", \"numberStacks\": 1, \"shape\": {\"aliases\": [], \"dimensions\": "
            "{\"A\": 0.032, \"B\": 0.0161, \"C\": 0.01065, \"D\": 0.0115, \"E\": 0.0232, \"F\": 0.0092, \"G\": 0.0, "
            "\"H\": 0.0}, \"family\": \"e\", \"familySubtype\": null, \"magneticCircuit\": \"open\", \"name\": \"E "
            "32/16/11\", \"type\": \"standard\"}, \"type\": \"two-piece set\"}, \"manufacturerInfo\": null, \"name\": "
            "\"My Core\"}"));
        auto winding = OpenMagnetics::CoilWrapper(
            json::parse("{\"bobbin\": \"Dummy\", \"functionalDescription\": [{\"isolationSide\": \"primary\", \"name\": "
                        "\"Primary\", \"numberParallels\": 1, \"numberTurns\": 33, \"wire\": \"Dummy\"}], "
                        "\"layersDescription\": null, \"sectionsDescription\": null, \"turnsDescription\": null}"));
        auto operatingPoint = OpenMagnetics::OperatingPoint(
            json::parse("{\"conditions\": {\"ambientRelativeHumidity\": null, \"ambientTemperature\": 37.0, "
                        "\"cooling\": null, \"name\": null}, \"excitationsPerWinding\": [{\"current\": {\"harmonics\": "
                        "null, \"processed\": null, \"waveform\": {\"ancillaryLabel\": null, \"data\": [-5.0, 5.0, "
                        "-5.0], \"numberPeriods\": null, \"time\": [0.0, 2.4999999999999998e-06, 1e-05]}}, "
                        "\"frequency\": 100000.0, \"magneticFieldStrength\": null, \"magneticFluxDensity\": null, "
                        "\"magnetizingCurrent\": null, \"name\": \"My Operating Point\"}], \"name\": null}"));

        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::string{models["gapReluctance"]});

        OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];

        auto magneticFluxDensity =
            magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint).second;

        excitation.set_magnetic_flux_density(magneticFluxDensity);
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();

        auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
        auto coreLosses = coreLossesModel->get_core_losses(core, excitation, temperature);

        CHECK_CLOSE(magneticFluxDensity.get_processed().value().get_offset(), 0, 0.0001);
    }
    TEST(Only_Current_1) {
        auto models = json::parse(R"({"coreLosses": "IGSE", "coreTemperature": "MANIKTALA", "gapReluctance": "ZHANG"})");
        auto core = OpenMagnetics::CoreWrapper(json::parse(R"({"functionalDescription": {"gapping": [{"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 0.001, "sectionDimensions": null, "shape": null, "type": "subtractive"}, {"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}, {"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}], "material": "3C97", "numberStacks": 1, "shape": {"aliases": [], "dimensions": {"A": 0.0391, "B": 0.0198, "C": 0.0125, "D": 0.0146, "E": 0.030100000000000002, "F": 0.0125, "G": 0.0, "H": 0.0}, "family": "etd", "familySubtype": "1", "magneticCircuit": null, "name": "ETD 39/20/13", "type": "standard"}, "type": "two-piece set"}, "manufacturerInfo": null, "name": "My Core"})"));
        auto winding = OpenMagnetics::CoilWrapper(
            json::parse(R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary", "numberParallels": 1, "numberTurns": 41, "wire": "Dummy"}], "layersDescription": null, "sectionsDescription": null, "turnsDescription": null})"));
        auto inputs = OpenMagnetics::InputsWrapper(
            json::parse(R"({"designRequirements": {"altitude": null, "cti": null, "insulationType": null, "leakageInductance": null, "magnetizingInductance": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.00034861070852064337}, "name": null, "operatingTemperature": null, "overvoltageCategory": null, "pollutionDegree": null, "turnsRatios": []}, "operatingPoints": [{"conditions": {"ambientRelativeHumidity": null, "ambientTemperature": 25.0, "cooling": null, "name": null}, "excitationsPerWinding": [{"current": {"harmonics": null, "processed": null, "waveform": {"ancillaryLabel": null, "data": [-8.0, 8.0, -8.0], "numberPeriods": null, "time": [0.0, 2.4999999999999998e-06, 1e-05]}}, "frequency": 100000.0, "magneticFieldStrength": null, "magneticFluxDensity": null, "magnetizingCurrent": null, "name": "My Operating Point", "voltage": null}], "name": null}]})"));

        auto operatingPoint = inputs.get_operating_point(0);
        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::string{models["gapReluctance"]});

        OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];

        auto magneticFluxDensity =
            magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint).second;

        excitation.set_magnetic_flux_density(magneticFluxDensity);
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();

        auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
        auto coreLosses = coreLossesModel->get_core_losses(core, excitation, temperature);

        CHECK_CLOSE(magneticFluxDensity.get_processed().value().get_offset(), 0, 0.0001);
    }
    TEST(Only_Voltage) {
        auto models =
            json::parse(R"({"coreLosses": "ROSHEN", "coreTemperature": "MANIKTALA", "gapReluctance": "ZHANG"})");
        auto core = OpenMagnetics::CoreWrapper(json::parse(
            R"({"functionalDescription": {"gapping": [{"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 0.001, "sectionDimensions": null, "shape": null, "type": "subtractive"}, {"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}, {"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}], "material": "3C97", "numberStacks": 1, "shape": {"aliases": [], "dimensions": {"A": 0.0391, "B": 0.0198, "C": 0.0125, "D": 0.0146, "E": 0.030100000000000002, "F": 0.0125, "G": 0.0, "H": 0.0}, "family": "etd", "familySubtype": "1", "magneticCircuit": null, "name": "ETD 39/20/13", "type": "standard"}, "type": "two-piece set"}, "manufacturerInfo": null, "name": "My Core"})"));
        auto winding = OpenMagnetics::CoilWrapper(json::parse(
            R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary", "numberParallels": 1, "numberTurns": 43, "wire": "Dummy"}], "layersDescription": null, "sectionsDescription": null, "turnsDescription": null})"));
        auto operatingPoint = OpenMagnetics::OperatingPoint(json::parse(
            R"({"conditions": {"ambientRelativeHumidity": null, "ambientTemperature": 25.0, "cooling": null, "name": null}, "excitationsPerWinding": [{"frequency": 100000.0, "magneticFieldStrength": null, "magneticFluxDensity": null, "magnetizingCurrent": null, "name": "My Operating Point", "voltage": {"harmonics": null, "processed": null, "waveform": {"ancillaryLabel": null, "data": [6.885, 6.885, -2.2949999999999995, -2.2949999999999995, 6.885], "numberPeriods": null, "time": [0.0, 2.4999999999999998e-06, 2.4999999999999998e-06, 1e-05, 1e-05]}}}], "name": null})"));

        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::string{models["gapReluctance"]});

        OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];

        auto magneticFluxDensity =
            magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint).second;

        excitation.set_magnetic_flux_density(magneticFluxDensity);
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();

        auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
        auto coreLosses = coreLossesModel->get_core_losses(core, excitation, temperature);

        CHECK_CLOSE(magneticFluxDensity.get_processed().value().get_offset(), 0, 0.0001);
    }
    TEST(Only_Current_With_Dc) {
        auto models = json::parse("{\"coreLosses\": \"IGSE\", \"gapReluctance\": \"BALAKRISHNAN\"}");
        auto core = OpenMagnetics::CoreWrapper(json::parse(
            "{\"functionalDescription\": {\"gapping\": [{\"area\": 9.8e-05, \"coordinates\": [0.0, 0.0001, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011301, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 0.0002, \"sectionDimensions\": [0.0092, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"subtractive\"}, {\"area\": 4.7e-05, \"coordinates\": [0.0138, 0.0, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011498, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.004401, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"residual\"}, {\"area\": 4.7e-05, \"coordinates\": [-0.0138, 0.0, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011498, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.004401, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"residual\"}], \"material\": \"3C97\", \"numberStacks\": 1, \"shape\": {\"aliases\": [], \"dimensions\": "
            "{\"A\": 0.032, \"B\": 0.0161, \"C\": 0.01065, \"D\": 0.0115, \"E\": 0.0232, \"F\": 0.0092, \"G\": 0.0, "
            "\"H\": 0.0}, \"family\": \"e\", \"familySubtype\": null, \"magneticCircuit\": \"open\", \"name\": \"E "
            "32/16/11\", \"type\": \"standard\"}, \"type\": \"two-piece set\"}, \"manufacturerInfo\": null, \"name\": "
            "\"My Core\"}"));
        auto winding = OpenMagnetics::CoilWrapper(
            json::parse("{\"bobbin\": \"Dummy\", \"functionalDescription\": [{\"isolationSide\": \"primary\", \"name\": "
                        "\"Primary\", \"numberParallels\": 1, \"numberTurns\": 33, \"wire\": \"Dummy\"}], "
                        "\"layersDescription\": null, \"sectionsDescription\": null, \"turnsDescription\": null}"));
        auto operatingPoint = OpenMagnetics::OperatingPoint(
            json::parse("{\"conditions\": {\"ambientRelativeHumidity\": null, \"ambientTemperature\": 37.0, "
                        "\"cooling\": null, \"name\": null}, \"excitationsPerWinding\": [{\"current\": {\"harmonics\": "
                        "null, \"processed\": null, \"waveform\": {\"ancillaryLabel\": null, \"data\": [-5.0, 15.0, "
                        "-5.0], \"numberPeriods\": null, \"time\": [0.0, 2.4999999999999998e-06, 1e-05]}}, "
                        "\"frequency\": 100000.0, \"magneticFieldStrength\": null, \"magneticFluxDensity\": null, "
                        "\"magnetizingCurrent\": null, \"name\": \"My Operating Point\"}], \"name\": null}"));

        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::string{models["gapReluctance"]});

        OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];

        auto magneticFluxDensity =
            magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint).second;

        excitation.set_magnetic_flux_density(magneticFluxDensity);
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();

        auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
        auto coreLosses = coreLossesModel->get_core_losses(core, excitation, temperature);

        CHECK_CLOSE(magneticFluxDensity.get_processed().value().get_offset(), 1, 1 * 0.1);
    }

    TEST(Crash_Voltage_Sin) {
        auto models = json::parse(R"({"coreLosses": "PROPRIETARY", "coreTemperature": "MANIKTALA", "gapReluctance": "ZHANG"})");
        auto core = OpenMagnetics::CoreWrapper(json::parse(R"({"functionalDescription": {"gapping": [{"length": 0.001, "type": "subtractive", "area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "sectionDimensions": null, "shape": null}, {"length": 1e-05, "type": "residual", "area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "sectionDimensions": null, "shape": null}, {"length": 1e-05, "type": "residual", "area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "sectionDimensions": null, "shape": null}], "material": "75-Series 60", "shape": {"family": "etd", "type": "standard", "aliases": [], "dimensions": {"A": 0.0391, "B": 0.0198, "C": 0.0125, "D": 0.0146, "E": 0.030100000000000002, "F": 0.0125, "G": 0.0, "H": 0.0}, "familySubtype": "1", "magneticCircuit": null, "name": "ETD 39/20/13"}, "type": "two-piece set", "numberStacks": 1}, "geometricalDescription": null, "manufacturerInfo": null, "name": "My Core", "processedDescription": null}      )"));
        auto winding = OpenMagnetics::CoilWrapper(
            json::parse(R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary", "numberParallels": 1, "numberTurns": 23, "wire": "Dummy"}], "layersDescription": null, "sectionsDescription": null, "turnsDescription": null})"));
        auto inputs = OpenMagnetics::InputsWrapper(
            json::parse(R"({"designRequirements": {"magnetizingInductance": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 3.8810245456770865e-05}, "turnsRatios": [], "altitude": null, "cti": null, "insulationType": null, "leakageInductance": null, "name": null, "operatingTemperature": null, "overvoltageCategory": null, "pollutionDegree": null}, "operatingPoints": [{"conditions": {"ambientTemperature": 25.0, "ambientRelativeHumidity": null, "cooling": null, "name": null}, "excitationsPerWinding": [{"frequency": 123000.0, "current": null, "magneticFieldStrength": null, "magneticFluxDensity": null, "magnetizingCurrent": null, "name": "My Operating Point Alf sin", "voltage": {"harmonics": null, "processed": null, "waveform": {"data": [0.0, 0.0, 15.4, 15.4, 0.0, 0.0, -15.4, -15.4, 0.0, 0.0], "numberPeriods": null, "ancillaryLabel": "Bipolar Rectangular", "time": [0.0, 1.2601626016260166e-06, 1.2601626016260166e-06, 2.8048780487804875e-06, 2.8048780487804875e-06, 5.325203252032522e-06, 5.325203252032522e-06, 6.869918699186992e-06, 6.869918699186992e-06, 8.130081300813007e-06]}}}], "name": null}]})"));

        auto operatingPoint = inputs.get_operating_point(0);
        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::string{models["gapReluctance"]});

        OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];

        auto magneticFluxDensity =
            magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint).second;

        excitation.set_magnetic_flux_density(magneticFluxDensity);
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();

        auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
        auto coreLosses = coreLossesModel->get_core_losses(core, excitation, temperature);

        CHECK_CLOSE(magneticFluxDensity.get_processed().value().get_offset(), 0, 0.001);
    }

    TEST(Crash_Losses) {
        auto models = json::parse(R"({"coreLosses": "STEINMETZ", "coreTemperature": "MANIKTALA", "gapReluctance": "ZHANG"})");
        auto core = OpenMagnetics::CoreWrapper(json::parse(R"({"functionalDescription": {"gapping": [{"length": 0.001, "type": "subtractive", "area": 0.000123, "coordinates": [0.0, 0.0005, 0.0], "distanceClosestNormalSurface": 0.0136, "distanceClosestParallelSurface": 0.0088, "sectionDimensions": [0.0125, 0.0125], "shape": "round"}, {"length": 5e-06, "type": "residual", "area": 6.2e-05, "coordinates": [0.017301, 0.0, 0.0], "distanceClosestNormalSurface": 0.014598, "distanceClosestParallelSurface": 0.0088, "sectionDimensions": [0.004501, 0.0125], "shape": "irregular"}, {"length": 5e-06, "type": "residual", "area": 6.2e-05, "coordinates": [-0.017301, 0.0, 0.0], "distanceClosestNormalSurface": 0.014598, "distanceClosestParallelSurface": 0.0088, "sectionDimensions": [0.004501, 0.0125], "shape": "irregular"}], "material": "3C95", "shape": {"family": "etd", "type": "standard", "aliases": [], "dimensions": {"A": 0.0391, "B": 0.0198, "C": 0.0125, "D": 0.0146, "E": 0.030100000000000002, "F": 0.0125, "G": 0.0, "H": 0.0}, "familySubtype": "1", "magneticCircuit": null, "name": "ETD 39/20/13"}, "type": "two-piece set", "numberStacks": 1}, "manufacturerInfo": null, "name": "My Core"})"));
        auto winding = OpenMagnetics::CoilWrapper(
            json::parse(R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary", "numberParallels": 1, "numberTurns": 10, "wire": "Dummy"}], "layersDescription": null, "sectionsDescription": null, "turnsDescription": null})"));
        auto inputs = OpenMagnetics::InputsWrapper(
            json::parse(R"({"designRequirements": {"magnetizingInductance": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 2e-05}, "turnsRatios": [], "altitude": null, "cti": null, "insulationType": null, "leakageInductance": null, "name": null, "operatingTemperature": null, "overvoltageCategory": null, "pollutionDegree": null}, "operatingPoints": [{"conditions": {"ambientTemperature": 25.0, "ambientRelativeHumidity": null, "cooling": null, "name": null}, "excitationsPerWinding": [{"frequency": 100000.0, "current": {"harmonics": null, "processed": null, "waveform": {"data": [-5.0, 5.0, -5.0], "numberPeriods": null, "ancillaryLabel": null, "time": [0.0, 2.5e-06, 1e-05]}}, "magneticFieldStrength": null, "magneticFluxDensity": null, "magnetizingCurrent": null, "name": "My Operating Point", "voltage": null}], "name": null}]})"));

        auto operatingPoint = inputs.get_operating_point(0);
        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::string{models["gapReluctance"]});

        OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];

        auto magneticFluxDensity =
            magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint).second;

        excitation.set_magnetic_flux_density(magneticFluxDensity);
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();

        auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
        auto coreLosses = coreLossesModel->get_core_losses(core, excitation, temperature);
    }

    TEST(Crash_Toroids) {
        auto models = json::parse(R"({"coreLosses": "IGSE", "coreTemperature": "MANIKTALA", "gapReluctance": "ZHANG"})");
        auto core = OpenMagnetics::CoreWrapper(json::parse(R"({"functionalDescription": {"gapping": [], "material": "3C97", "shape": {"family": "t", "type": "standard", "aliases": [], "dimensions": {"A": 0.0034300000000000003, "B": 0.0017800000000000001, "C": 0.00203}, "familySubtype": null, "magneticCircuit": "closed", "name": "T 3.43/1.78/2.03"}, "type": "toroidal", "numberStacks": 1}, "geometricalDescription": null, "manufacturerInfo": null, "name": "My Core", "processedDescription": null})"));
        auto winding = OpenMagnetics::CoilWrapper(
            json::parse(R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary", "numberParallels": 1, "numberTurns": 26, "wire": "Dummy"}], "layersDescription": null, "sectionsDescription": null, "turnsDescription": null})"));
        auto inputs = OpenMagnetics::InputsWrapper(
            json::parse(R"({"designRequirements": {"magnetizingInductance": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.0009443757859214556}, "turnsRatios": [], "altitude": null, "cti": null, "insulationType": null, "leakageInductance": null, "name": null, "operatingTemperature": null, "overvoltageCategory": null, "pollutionDegree": null}, "operatingPoints": [{"conditions": {"ambientTemperature": 25.0, "ambientRelativeHumidity": null, "cooling": null, "name": null}, "excitationsPerWinding": [{"frequency": 100000.0, "current": {"harmonics": null, "processed": null, "waveform": {"data": [-5.0, 5.0, -5.0], "numberPeriods": null, "ancillaryLabel": null, "time": [0.0, 2.5e-06, 1e-05]}}, "magneticFieldStrength": null, "magneticFluxDensity": null, "magnetizingCurrent": null, "name": "My Operating Point", "voltage": null}], "name": null}]})"));

        auto operatingPoint = inputs.get_operating_point(0);
        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::string{models["gapReluctance"]});

        OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];

        auto magneticFluxDensity =
            magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint).second;

        excitation.set_magnetic_flux_density(magneticFluxDensity);
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();

        auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
        auto coreLosses = coreLossesModel->get_core_losses(core, excitation, temperature);

        CHECK_CLOSE(magneticFluxDensity.get_processed().value().get_offset(), 0, 0.0001);
    }

    TEST(Test_Methods) {
        std::vector<std::string> methods = OpenMagnetics::CoreLossesModel::get_methods_string("3C97");
        CHECK(methods.size() == 6);
        CHECK(methods[0] == "steinmetz");
        CHECK(methods[1] == "igse");
        CHECK(methods[2] == "barg");
        CHECK(methods[3] == "albach");
        CHECK(methods[4] == "mse");
        CHECK(methods[5] == "roshen");
        methods = OpenMagnetics::CoreLossesModel::get_methods_string("XFlux 19");
        CHECK(methods.size() == 1);
        CHECK(methods[0] == "proprietary");
        methods = OpenMagnetics::CoreLossesModel::get_methods_string("A07");
        CHECK(methods.size() == 1);
        CHECK(methods[0] == "loss_factor");
    }

    TEST(Test_Manufacturer_Magnetics) {
        auto models =
            json::parse(R"({"coreLosses": "PROPRIETARY", "coreTemperature": "MANIKTALA", "gapReluctance": "ZHANG"})");
        auto core = OpenMagnetics::CoreWrapper(json::parse(
            R"({"functionalDescription": {"gapping": [{"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 0.001, "sectionDimensions": null, "shape": null, "type": "subtractive"}, {"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}, {"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}], "material": "XFlux 60", "numberStacks": 1, "shape": {"aliases": [], "dimensions": {"A": 0.0391, "B": 0.0198, "C": 0.0125, "D": 0.0146, "E": 0.030100000000000002, "F": 0.0125, "G": 0.0, "H": 0.0}, "family": "etd", "familySubtype": "1", "magneticCircuit": null, "name": "ETD 39/20/13", "type": "standard"}, "type": "two-piece set"}, "manufacturerInfo": null, "name": "My Core"})"));
        auto winding = OpenMagnetics::CoilWrapper(json::parse(
            R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary", "numberParallels": 1, "numberTurns": 43, "wire": "Dummy"}], "layersDescription": null, "sectionsDescription": null, "turnsDescription": null})"));
        auto operatingPoint = OpenMagnetics::OperatingPoint(json::parse(
            R"({"conditions": {"ambientRelativeHumidity": null, "ambientTemperature": 25.0, "cooling": null, "name": null}, "excitationsPerWinding": [{"frequency": 100000.0, "magneticFieldStrength": null, "magneticFluxDensity": null, "magnetizingCurrent": null, "name": "My Operating Point", "voltage": {"harmonics": null, "processed": null, "waveform": {"ancillaryLabel": null, "data": [688.5, 688.5, -229.49999999999995, -229.49999999999995, 688.5], "numberPeriods": null, "time": [0.0, 2.4999999999999998e-06, 2.4999999999999998e-06, 1e-05, 1e-05]}}}], "name": null})"));

        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::string{models["gapReluctance"]});

        OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];

        auto magneticFluxDensity =
            magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint).second;

        excitation.set_magnetic_flux_density(magneticFluxDensity);
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();

        auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
        auto coreLosses = coreLossesModel->get_core_losses(core, excitation, temperature);

        CHECK_CLOSE(magneticFluxDensity.get_processed().value().get_offset(), 0, 0.002);
        CHECK_CLOSE(coreLosses.get_core_losses(), 31.2, 31.2 * maxError);
    }

    TEST(Test_Manufacturer_Micrometals) {
        auto models =
            json::parse(R"({"coreLosses": "PROPRIETARY", "coreTemperature": "MANIKTALA", "gapReluctance": "ZHANG"})");
        auto core = OpenMagnetics::CoreWrapper(json::parse(
            R"({"functionalDescription": {"gapping": [{"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 0.001, "sectionDimensions": null, "shape": null, "type": "subtractive"}, {"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}, {"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}], "material": "MS 40", "numberStacks": 1, "shape": {"aliases": [], "dimensions": {"A": 0.0391, "B": 0.0198, "C": 0.0125, "D": 0.0146, "E": 0.030100000000000002, "F": 0.0125, "G": 0.0, "H": 0.0}, "family": "etd", "familySubtype": "1", "magneticCircuit": null, "name": "ETD 39/20/13", "type": "standard"}, "type": "two-piece set"}, "manufacturerInfo": null, "name": "My Core"})"));
        auto winding = OpenMagnetics::CoilWrapper(json::parse(
            R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary", "numberParallels": 1, "numberTurns": 43, "wire": "Dummy"}], "layersDescription": null, "sectionsDescription": null, "turnsDescription": null})"));
        auto operatingPoint = OpenMagnetics::OperatingPoint(json::parse(
            R"({"conditions": {"ambientRelativeHumidity": null, "ambientTemperature": 25.0, "cooling": null, "name": null}, "excitationsPerWinding": [{"frequency": 100000.0, "magneticFieldStrength": null, "magneticFluxDensity": null, "magnetizingCurrent": null, "name": "My Operating Point", "voltage": {"harmonics": null, "processed": null, "waveform": {"ancillaryLabel": null, "data": [688.5, 688.5, -229.49999999999995, -229.49999999999995, 688.5], "numberPeriods": null, "time": [0.0, 2.4999999999999998e-06, 2.4999999999999998e-06, 1e-05, 1e-05]}}}], "name": null})"));

        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::string{models["gapReluctance"]});

        OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];

        auto magneticFluxDensity =
            magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint).second;

        excitation.set_magnetic_flux_density(magneticFluxDensity);
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();

        auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
        auto coreLosses = coreLossesModel->get_core_losses(core, excitation, temperature);

        CHECK_CLOSE(magneticFluxDensity.get_processed().value().get_offset(), 0, 0.002);
        CHECK_CLOSE(coreLosses.get_core_losses(), 23.1, 23.1 * maxError);
    }

    TEST(Test_XFlux_19) {
        auto models = json::parse("{\"coreLosses\": \"PROPRIETARY\", \"gapReluctance\": \"BALAKRISHNAN\"}");
        auto core = OpenMagnetics::CoreWrapper(json::parse(
            "{\"functionalDescription\": {\"gapping\": [{\"area\": 9.8e-05, \"coordinates\": [0.0, 0.0001, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011301, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 0.0002, \"sectionDimensions\": [0.0092, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"subtractive\"}, {\"area\": 4.7e-05, \"coordinates\": [0.0138, 0.0, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011498, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.004401, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"residual\"}, {\"area\": 4.7e-05, \"coordinates\": [-0.0138, 0.0, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011498, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.004401, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"residual\"}], \"material\": \"XFlux 19\", \"numberStacks\": 1, \"shape\": {\"aliases\": [], "
            "\"dimensions\": {\"A\": 0.032, \"B\": 0.0161, \"C\": 0.01065, \"D\": 0.0115, \"E\": 0.0232, \"F\": "
            "0.0092, \"G\": 0.0, \"H\": 0.0}, \"family\": \"e\", \"familySubtype\": null, \"magneticCircuit\": "
            "\"open\", \"name\": \"E 32/16/11\", \"type\": \"standard\"}, \"type\": \"two-piece set\"}, "
            "\"manufacturerInfo\": null, \"name\": \"My Core\"}"));
        auto winding = OpenMagnetics::CoilWrapper(
            json::parse("{\"bobbin\": \"Dummy\", \"functionalDescription\": [{\"isolationSide\": \"primary\", \"name\": "
                        "\"Primary\", \"numberParallels\": 1, \"numberTurns\": 33, \"wire\": \"Dummy\"}], "
                        "\"layersDescription\": null, \"sectionsDescription\": null, \"turnsDescription\": null}"));
        auto operatingPoint = OpenMagnetics::OperatingPoint(
            json::parse("{\"conditions\": {\"ambientRelativeHumidity\": null, \"ambientTemperature\": 37.0, "
                        "\"cooling\": null, \"name\": null}, \"excitationsPerWinding\": [{\"current\": {\"harmonics\": "
                        "null, \"processed\": null, \"waveform\": {\"ancillaryLabel\": null, \"data\": [-5.0, 5.0, "
                        "-5.0], \"numberPeriods\": null, \"time\": [0.0, 2.4999999999999998e-06, 1e-05]}}, "
                        "\"frequency\": 100000.0, \"magneticFieldStrength\": null, \"magneticFluxDensity\": null, "
                        "\"magnetizingCurrent\": null, \"name\": \"My Operating Point\"}], \"name\": null}"));

        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::string{models["gapReluctance"]});

        OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];

        auto magneticFluxDensity =
            magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint).second;

        excitation.set_magnetic_flux_density(magneticFluxDensity);
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();

        auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
        auto coreLosses = coreLossesModel->get_core_losses(core, excitation, temperature);
        auto calculatedVolumetricCoreLosses = coreLosses.get_volumetric_losses().value();
        auto expectedVolumetricLosses = 297420;

        double maximumAdmittedErrorVolumetricCoreLossesValue =
            maximumAdmittedErrorVolumetricCoreLosses[OpenMagnetics::CoreLossesModels::STEINMETZ];
        CHECK_CLOSE(magneticFluxDensity.get_processed().value().get_offset(), 0, 0.0001);
        CHECK_CLOSE(calculatedVolumetricCoreLosses, expectedVolumetricLosses,
                    expectedVolumetricLosses * maximumAdmittedErrorVolumetricCoreLossesValue);
    }

    TEST(Test_A07) {
        auto models = json::parse("{\"coreLosses\": \"LOSS_FACTOR\", \"gapReluctance\": \"BALAKRISHNAN\"}");
        auto core = OpenMagnetics::CoreWrapper(json::parse(R"({"distributorsInfo": [], "functionalDescription": {"coating": null, "gapping": [{"area": 0.000078, "coordinates": [0, 0.0003, 0 ], "distanceClosestNormalSurface": 0.00865, "distanceClosestParallelSurface": 0.005325, "length": 0.0006, "sectionDimensions": [0.00725, 0.01075 ], "shape": "rectangular", "type": "subtractive" }, {"area": 0.000039, "coordinates": [0.010738, 0, 0 ], "distanceClosestNormalSurface": 0.00895, "distanceClosestParallelSurface": 0.005325, "length": 0.000005, "sectionDimensions": [0.003575, 0.01075 ], "shape": "rectangular", "type": "residual" }, {"area": 0.000039, "coordinates": [-0.010738, 0, 0 ], "distanceClosestNormalSurface": 0.00895, "distanceClosestParallelSurface": 0.005325, "length": 0.000005, "sectionDimensions": [0.003575, 0.01075 ], "shape": "rectangular", "type": "residual" } ], "material": "A07", "numberStacks": 1, "shape": {"aliases": [], "dimensions": {"A": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0258, "minimum": 0.0243, "nominal": null }, "B": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0128, "minimum": 0.0123, "nominal": null }, "C": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.011, "minimum": 0.0105, "nominal": null }, "D": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0092, "minimum": 0.0087, "nominal": null }, "E": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0183, "minimum": 0.0175, "nominal": null }, "F": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0075, "minimum": 0.007, "nominal": null } }, "family": "e", "familySubtype": null, "magneticCircuit": "open", "name": "E 25/13/11", "type": "standard" }, "type": "two-piece set", "magneticCircuit": "open" }, "geometricalDescription": [{"coordinates": [0, 0, 0 ], "dimensions": null, "insulationMaterial": null, "machining": [{"coordinates": [0, 0.0003, 0 ], "length": 0.0006 } ], "material": "A07", "rotation": [3.141592653589793, 3.141592653589793, 0 ], "shape": {"aliases": [], "dimensions": {"A": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0258, "minimum": 0.0243, "nominal": null }, "B": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0128, "minimum": 0.0123, "nominal": null }, "C": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.011, "minimum": 0.0105, "nominal": null }, "D": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0092, "minimum": 0.0087, "nominal": null }, "E": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0183, "minimum": 0.0175, "nominal": null }, "F": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0075, "minimum": 0.007, "nominal": null } }, "family": "e", "familySubtype": null, "magneticCircuit": "open", "name": "E 25/13/11", "type": "standard" }, "type": "half set" }, {"coordinates": [0, 0, 0 ], "dimensions": null, "insulationMaterial": null, "machining": null, "material": "A07", "rotation": [0, 0, 0 ], "shape": {"aliases": [], "dimensions": {"A": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0258, "minimum": 0.0243, "nominal": null }, "B": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0128, "minimum": 0.0123, "nominal": null }, "C": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.011, "minimum": 0.0105, "nominal": null }, "D": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0092, "minimum": 0.0087, "nominal": null }, "E": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0183, "minimum": 0.0175, "nominal": null }, "F": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0075, "minimum": 0.007, "nominal": null } }, "family": "e", "familySubtype": null, "magneticCircuit": "open", "name": "E 25/13/11", "type": "standard" }, "type": "half set" } ], "manufacturerInfo": null, "name": "Custom", "processedDescription": {"columns": [{"area": 0.000078, "coordinates": [0, 0, 0 ], "depth": 0.01075, "height": 0.0179, "minimumDepth": null, "minimumWidth": null, "shape": "rectangular", "type": "central", "width": 0.00725 }, {"area": 0.000039, "coordinates": [0.010738, 0, 0 ], "depth": 0.01075, "height": 0.0179, "minimumDepth": null, "minimumWidth": null, "shape": "rectangular", "type": "lateral", "width": 0.003575 }, {"area": 0.000039, "coordinates": [-0.010738, 0, 0 ], "depth": 0.01075, "height": 0.0179, "minimumDepth": null, "minimumWidth": null, "shape": "rectangular", "type": "lateral", "width": 0.003575 } ], "depth": 0.01075, "effectiveParameters": {"effectiveArea": 0.00007739519022938956, "effectiveLength": 0.05775787070464925, "effectiveVolume": 0.000004470181390430815, "minimumArea": 0.00007686249999999999 }, "height": 0.0251, "width": 0.02505, "windingWindows": [{"angle": null, "area": 0.00009531749999999999, "coordinates": [0.0036249999999999998, 0 ], "height": 0.0179, "radialHeight": null, "sectionsAlignment": null, "sectionsOrientation": null, "shape": null, "width": 0.005325 } ] } })"));
        auto winding = OpenMagnetics::CoilWrapper(json::parse(R"({"bobbin": {"distributorsInfo": null, "functionalDescription": null, "manufacturerInfo": null, "name": null, "processedDescription": {"columnDepth": 0.0065, "columnShape": "rectangular", "columnThickness": 0.0011250000000000001, "columnWidth": 0.004750000000000001, "coordinates": [0, 0, 0 ], "pins": null, "wallThickness": 0.0010499999999999989, "windingWindows": [{"angle": null, "area": 0.00006636000000000001, "coordinates": [0.00685, 0, 0 ], "height": 0.0158, "radialHeight": null, "sectionsAlignment": "inner or top", "sectionsOrientation": "overlapping", "shape": "rectangular", "width": 0.0042 } ] } }, "functionalDescription": [{"connections": null, "isolationSide": "primary", "name": "Primary", "numberParallels": 1, "numberTurns": 32, "wire": {"coating": {"breakdownVoltage": 2700, "grade": 1, "material": null, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "enamelled" }, "conductingArea": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 3.6637960384511227e-7 }, "conductingDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.000683 }, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Nearson", "orderCode": null, "reference": null, "status": null }, "material": "copper", "name": "Round 21.5 - Single Build", "numberConductors": 1, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.000716 }, "outerHeight": null, "outerWidth": null, "standard": "NEMA MW 1000 C", "standardName": "21.5 AWG", "strand": null, "type": "round" } } ], "layersDescription": [{"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.0051080000000000006, 0 ], "dimensions": [0.000716, 0.015528500000000004 ], "fillingFactor": 0.7250632911392404, "insulationMaterial": null, "name": "Primary section 0 layer 0", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [0.5 ], "winding": "Primary" } ], "section": "Primary section 0", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveTurns" }, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.005824000000000001, 0 ], "dimensions": [0.000716, 0.015528500000000004 ], "fillingFactor": 0.7250632911392404, "insulationMaterial": null, "name": "Primary section 0 layer 1", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [0.5 ], "winding": "Primary" } ], "section": "Primary section 0", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveTurns" } ], "sectionsDescription": [{"coordinateSystem": "cartesian", "coordinates": [0.005466000000000001, 0 ], "dimensions": [0.0014320000000000003, 0.015528500000000004 ], "fillingFactor": 0.24721205545509337, "layersAlignment": null, "layersOrientation": "overlapping", "margin": [0, 0 ], "name": "Primary section 0", "partialWindings": [{"connections": null, "parallelsProportion": [1 ], "winding": "Primary" } ], "type": "conduction", "windingStyle": "windByConsecutiveTurns" } ], "turnsDescription": [{"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051080000000000006, 0.007406250000000002 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 0", "length": 0.04724938033997029, "name": "Primary parallel 0 turn 0", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051080000000000006, 0.006418750000000003 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 0", "length": 0.04724938033997029, "name": "Primary parallel 0 turn 1", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051080000000000006, 0.005431250000000002 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 0", "length": 0.04724938033997029, "name": "Primary parallel 0 turn 2", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051080000000000006, 0.0044437500000000015 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 0", "length": 0.04724938033997029, "name": "Primary parallel 0 turn 3", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051080000000000006, 0.0034562500000000014 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 0", "length": 0.04724938033997029, "name": "Primary parallel 0 turn 4", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051080000000000006, 0.0024687500000000013 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 0", "length": 0.04724938033997029, "name": "Primary parallel 0 turn 5", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051080000000000006, 0.0014812500000000012 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 0", "length": 0.04724938033997029, "name": "Primary parallel 0 turn 6", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051080000000000006, 0.0004937500000000011 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 0", "length": 0.04724938033997029, "name": "Primary parallel 0 turn 7", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051080000000000006, -0.000493749999999999 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 0", "length": 0.04724938033997029, "name": "Primary parallel 0 turn 8", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051080000000000006, -0.001481249999999999 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 0", "length": 0.04724938033997029, "name": "Primary parallel 0 turn 9", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051080000000000006, -0.002468749999999999 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 0", "length": 0.04724938033997029, "name": "Primary parallel 0 turn 10", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051080000000000006, -0.0034562499999999993 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 0", "length": 0.04724938033997029, "name": "Primary parallel 0 turn 11", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051080000000000006, -0.00444375 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 0", "length": 0.04724938033997029, "name": "Primary parallel 0 turn 12", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051080000000000006, -0.00543125 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 0", "length": 0.04724938033997029, "name": "Primary parallel 0 turn 13", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051080000000000006, -0.006418750000000001 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 0", "length": 0.04724938033997029, "name": "Primary parallel 0 turn 14", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051080000000000006, -0.007406250000000001 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 0", "length": 0.04724938033997029, "name": "Primary parallel 0 turn 15", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005824000000000001, 0.007406250000000002 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 1", "length": 0.051748141019910876, "name": "Primary parallel 0 turn 16", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005824000000000001, 0.006418750000000003 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 1", "length": 0.051748141019910876, "name": "Primary parallel 0 turn 17", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005824000000000001, 0.005431250000000002 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 1", "length": 0.051748141019910876, "name": "Primary parallel 0 turn 18", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005824000000000001, 0.0044437500000000015 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 1", "length": 0.051748141019910876, "name": "Primary parallel 0 turn 19", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005824000000000001, 0.0034562500000000014 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 1", "length": 0.051748141019910876, "name": "Primary parallel 0 turn 20", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005824000000000001, 0.0024687500000000013 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 1", "length": 0.051748141019910876, "name": "Primary parallel 0 turn 21", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005824000000000001, 0.0014812500000000012 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 1", "length": 0.051748141019910876, "name": "Primary parallel 0 turn 22", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005824000000000001, 0.0004937500000000011 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 1", "length": 0.051748141019910876, "name": "Primary parallel 0 turn 23", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005824000000000001, -0.000493749999999999 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 1", "length": 0.051748141019910876, "name": "Primary parallel 0 turn 24", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005824000000000001, -0.001481249999999999 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 1", "length": 0.051748141019910876, "name": "Primary parallel 0 turn 25", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005824000000000001, -0.002468749999999999 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 1", "length": 0.051748141019910876, "name": "Primary parallel 0 turn 26", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005824000000000001, -0.0034562499999999993 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 1", "length": 0.051748141019910876, "name": "Primary parallel 0 turn 27", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005824000000000001, -0.00444375 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 1", "length": 0.051748141019910876, "name": "Primary parallel 0 turn 28", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005824000000000001, -0.00543125 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 1", "length": 0.051748141019910876, "name": "Primary parallel 0 turn 29", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005824000000000001, -0.006418750000000001 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 1", "length": 0.051748141019910876, "name": "Primary parallel 0 turn 30", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" }, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005824000000000001, -0.007406250000000001 ], "dimensions": [0.000716, 0.000716 ], "layer": "Primary section 0 layer 1", "length": 0.051748141019910876, "name": "Primary parallel 0 turn 31", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary" } ] })"));
        auto operatingPoint = OpenMagnetics::OperatingPoint(json::parse(R"({"name": "Operating Point No. 1", "conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5 ], "time": [0, 0.000005, 0.00001 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular", "acEffectiveFrequency": 110746.40291779551, "effectiveFrequency": 110746.40291779551, "peak": 5, "rms": 2.8874560332150576, "thd": 0.12151487440704967 }, "harmonics": {"amplitudes": [1.1608769501236793e-14, 4.05366124583194, 1.787369544444173e-15, 0.4511310569983995, 9.749053004706756e-16, 0.16293015292554872, 4.036157626725542e-16, 0.08352979924600704, 3.4998295008010614e-16, 0.0508569581336163, 3.1489164048780735e-16, 0.034320410449418075, 3.142469873118059e-16, 0.024811988673843106, 2.3653352035940994e-16, 0.018849001010678823, 2.9306524147249266e-16, 0.014866633059596499, 1.796485796132569e-16, 0.012077180559557785, 1.6247782523152451e-16, 0.010049063750920609, 1.5324769149805092e-16, 0.008529750975091871, 1.0558579733068502e-16, 0.007363501410705499, 7.513269775674661e-17, 0.006450045785294609, 5.871414177162291e-17, 0.005722473794997712, 9.294731722001391e-17, 0.005134763398167541, 1.194820309200107e-16, 0.004654430423785411, 8.2422739080512e-17, 0.004258029771397032, 9.5067306351894e-17, 0.0039283108282380024, 1.7540347128474968e-16, 0.0036523670873925395, 9.623794010508822e-17, 0.0034204021424253787, 1.4083470894369491e-16, 0.0032248884817922927, 1.4749333016985644e-16, 0.0030599828465501895, 1.0448590642474364e-16, 0.002921112944200383, 7.575487373767413e-18, 0.002804680975178716, 7.419510610361002e-17, 0.0027078483284668376, 3.924741709073613e-17, 0.0026283777262804953, 2.2684279102637236e-17, 0.0025645167846443107, 8.997077625295079e-17, 0.0025149120164513483, 7.131074184849219e-17, 0.0024785457043284276, 9.354417496250849e-17, 0.0024546904085875065, 1.2488589642405877e-16, 0.0024428775264784264 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] } }, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5 ], "time": [0, 0, 0.000005, 0.000005, 0.00001 ] }, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 591485.5360118389, "effectiveFrequency": 553357.3374711702, "peak": 70.5, "rms": 51.572309672924284, "thd": 0.4833151484524849 }, "harmonics": {"amplitudes": [24.2890625, 57.92076613061847, 1.421875, 19.27588907896988, 1.421875, 11.528257939603122, 1.421875, 8.194467538528329, 1.421875, 6.331896912839248, 1.421875, 5.137996046859012, 1.421875, 4.304077056139349, 1.421875, 3.6860723299088454, 1.421875, 3.207698601961777, 1.421875, 2.8247804703632298, 1.421875, 2.509960393415185, 1.421875, 2.2453859950684323, 1.421875, 2.01890737840567, 1.421875, 1.8219644341144872, 1.421875, 1.6483482744897402, 1.421875, 1.4934420157473332, 1.421875, 1.3537375367153817, 1.421875, 1.2265178099275544, 1.421875, 1.1096421410704556, 1.421875, 1.0013973584174929, 1.421875, 0.9003924136274832, 1.421875, 0.8054822382572133, 1.421875, 0.7157117294021269, 1.421875, 0.6302738400635857, 1.421875, 0.5484777114167545, 1.421875, 0.46972405216147894, 1.421875, 0.3934858059809043, 1.421875, 0.31929270856030145, 1.421875, 0.24671871675852053, 1.421875, 0.17537155450693565, 1.421875, 0.10488380107099537, 1.421875, 0.034905072061178544 ], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000 ] } } } ]})"));

        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::string{models["gapReluctance"]});

        OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];

        auto magneticFluxDensity =
            magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint).second;

        excitation.set_magnetic_flux_density(magneticFluxDensity);
        excitation.set_magnetizing_current(excitation.get_current().value());
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();

        auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
        auto coreLosses = coreLossesModel->get_core_losses(core, excitation, temperature);
        auto calculatedVolumetricCoreLosses = coreLosses.get_volumetric_losses().value();
        auto expectedVolumetricLosses = 297420;

        double maximumAdmittedErrorVolumetricCoreLossesValue =
            maximumAdmittedErrorVolumetricCoreLosses[OpenMagnetics::CoreLossesModels::STEINMETZ];
        CHECK_CLOSE(magneticFluxDensity.get_processed().value().get_offset(), 0, 0.0001);
        CHECK_CLOSE(calculatedVolumetricCoreLosses, expectedVolumetricLosses,
                    expectedVolumetricLosses * maximumAdmittedErrorVolumetricCoreLossesValue);
    }
}

SUITE(FrequencyFromCoreLosses) {
    double maxError = 0.05;
    TEST(Frequency_From_Losses_Steinmetz) {
        auto models = json::parse("{\"coreLosses\": \"STEINMETZ\", \"gapReluctance\": \"BALAKRISHNAN\"}");
        auto core = OpenMagnetics::CoreWrapper(json::parse(
            "{\"functionalDescription\": {\"gapping\": [{\"area\": 9.8e-05, \"coordinates\": [0.0, 0.0001, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011301, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 0.0002, \"sectionDimensions\": [0.0092, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"subtractive\"}, {\"area\": 4.7e-05, \"coordinates\": [0.0138, 0.0, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011498, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.004401, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"residual\"}, {\"area\": 4.7e-05, \"coordinates\": [-0.0138, 0.0, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011498, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.004401, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"residual\"}], \"material\": \"3C97\", \"numberStacks\": 1, \"shape\": {\"aliases\": [], \"dimensions\": "
            "{\"A\": 0.032, \"B\": 0.0161, \"C\": 0.01065, \"D\": 0.0115, \"E\": 0.0232, \"F\": 0.0092, \"G\": 0.0, "
            "\"H\": 0.0}, \"family\": \"e\", \"familySubtype\": null, \"magneticCircuit\": \"open\", \"name\": \"E "
            "32/16/11\", \"type\": \"standard\"}, \"type\": \"two-piece set\"}, \"manufacturerInfo\": null, \"name\": "
            "\"My Core\"}"));
        auto winding = OpenMagnetics::CoilWrapper(
            json::parse("{\"bobbin\": \"Dummy\", \"functionalDescription\": [{\"isolationSide\": \"primary\", \"name\": "
                        "\"Primary\", \"numberParallels\": 1, \"numberTurns\": 33, \"wire\": \"Dummy\"}], "
                        "\"layersDescription\": null, \"sectionsDescription\": null, \"turnsDescription\": null}"));
        auto operatingPoint = OpenMagnetics::OperatingPoint(json::parse(
            "{\"conditions\": {\"ambientRelativeHumidity\": null, \"ambientTemperature\": 37.0, \"cooling\": null, "
            "\"name\": null}, \"excitationsPerWinding\": [{\"current\": null, \"frequency\": 100000.0, \"magneticFieldStrength\": null, "
            "\"magneticFluxDensity\": null, \"magnetizingCurrent\": null, \"name\": \"My Operating Point\", "
            "\"voltage\": {\"harmonics\": null, \"processed\": null, \"waveform\": {\"ancillaryLabel\": null, "
            "\"data\": [7.5, 7.5, -2.4999999999999996, -2.4999999999999996, 7.5], \"numberPeriods\": null, \"time\": "
            "[0.0, 2.4999999999999998e-06, 2.4999999999999998e-06, 1e-05, 1e-05]}}}], \"name\": null}"));

        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::string{models["gapReluctance"]});
        OpenMagnetics::InputsWrapper::scale_time_to_frequency(operatingPoint, 324578);

        OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];

        auto magneticFluxDensity =
            magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint).second;

        excitation.set_magnetic_flux_density(magneticFluxDensity);
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();


        auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
        auto coreLosses = coreLossesModel->get_core_losses(core, excitation, temperature);


        auto frequencyFromCoreLosses = coreLossesModel->get_frequency_from_core_losses(core, magneticFluxDensity, temperature, coreLosses.get_core_losses()); 
        CHECK_CLOSE(excitation.get_frequency(), frequencyFromCoreLosses, frequencyFromCoreLosses * maxError);

    }
    TEST(Frequency_From_Losses_Igse) {
        auto models = json::parse("{\"coreLosses\": \"IGSE\", \"gapReluctance\": \"BALAKRISHNAN\"}");
        auto core = OpenMagnetics::CoreWrapper(json::parse(
            "{\"functionalDescription\": {\"gapping\": [{\"area\": 9.8e-05, \"coordinates\": [0.0, 0.0001, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011301, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 0.0002, \"sectionDimensions\": [0.0092, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"subtractive\"}, {\"area\": 4.7e-05, \"coordinates\": [0.0138, 0.0, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011498, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.004401, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"residual\"}, {\"area\": 4.7e-05, \"coordinates\": [-0.0138, 0.0, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011498, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.004401, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"residual\"}], \"material\": \"3C97\", \"numberStacks\": 1, \"shape\": {\"aliases\": [], \"dimensions\": "
            "{\"A\": 0.032, \"B\": 0.0161, \"C\": 0.01065, \"D\": 0.0115, \"E\": 0.0232, \"F\": 0.0092, \"G\": 0.0, "
            "\"H\": 0.0}, \"family\": \"e\", \"familySubtype\": null, \"magneticCircuit\": \"open\", \"name\": \"E "
            "32/16/11\", \"type\": \"standard\"}, \"type\": \"two-piece set\"}, \"manufacturerInfo\": null, \"name\": "
            "\"My Core\"}"));
        auto winding = OpenMagnetics::CoilWrapper(
            json::parse("{\"bobbin\": \"Dummy\", \"functionalDescription\": [{\"isolationSide\": \"primary\", \"name\": "
                        "\"Primary\", \"numberParallels\": 1, \"numberTurns\": 33, \"wire\": \"Dummy\"}], "
                        "\"layersDescription\": null, \"sectionsDescription\": null, \"turnsDescription\": null}"));
        auto operatingPoint = OpenMagnetics::OperatingPoint(json::parse(
            "{\"conditions\": {\"ambientRelativeHumidity\": null, \"ambientTemperature\": 37.0, \"cooling\": null, "
            "\"name\": null}, \"excitationsPerWinding\": [{\"current\": {\"harmonics\": null, \"processed\": null, "
            "\"waveform\": {\"ancillaryLabel\": null, \"data\": [-5.0, 5.0, -5.0], \"numberPeriods\": null, \"time\": "
            "[0.0, 2.4999999999999998e-06, 1e-05]}}, \"frequency\": 100000.0, \"magneticFieldStrength\": null, "
            "\"magneticFluxDensity\": null, \"magnetizingCurrent\": null, \"name\": \"My Operating Point\", "
            "\"voltage\": {\"harmonics\": null, \"processed\": null, \"waveform\": {\"ancillaryLabel\": null, "
            "\"data\": [7.5, 7.5, -2.4999999999999996, -2.4999999999999996, 7.5], \"numberPeriods\": null, \"time\": "
            "[0.0, 2.4999999999999998e-06, 2.4999999999999998e-06, 1e-05, 1e-05]}}}], \"name\": null}"));

        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::string{models["gapReluctance"]});
        OpenMagnetics::InputsWrapper::scale_time_to_frequency(operatingPoint, 324578);

        OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];

        auto magneticFluxDensity =
            magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint).second;

        excitation.set_magnetic_flux_density(magneticFluxDensity);
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();

        auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
        auto coreLosses = coreLossesModel->get_core_losses(core, excitation, temperature);


        auto frequencyFromCoreLosses = coreLossesModel->get_frequency_from_core_losses(core, magneticFluxDensity, temperature, coreLosses.get_core_losses()); 
        CHECK_CLOSE(excitation.get_frequency(), frequencyFromCoreLosses, frequencyFromCoreLosses * maxError);

    }
    TEST(Frequency_From_Losses_Magnetics) {
        auto models = json::parse("{\"coreLosses\": \"PROPRIETARY\", \"gapReluctance\": \"BALAKRISHNAN\"}");
        auto core = OpenMagnetics::CoreWrapper(json::parse(
            "{\"functionalDescription\": {\"gapping\": [{\"area\": 9.8e-05, \"coordinates\": [0.0, 0.0001, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011301, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 0.0002, \"sectionDimensions\": [0.0092, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"subtractive\"}, {\"area\": 4.7e-05, \"coordinates\": [0.0138, 0.0, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011498, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.004401, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"residual\"}, {\"area\": 4.7e-05, \"coordinates\": [-0.0138, 0.0, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011498, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.004401, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"residual\"}], \"material\": \"XFlux 19\", \"numberStacks\": 1, \"shape\": {\"aliases\": [], "
            "\"dimensions\": {\"A\": 0.032, \"B\": 0.0161, \"C\": 0.01065, \"D\": 0.0115, \"E\": 0.0232, \"F\": "
            "0.0092, \"G\": 0.0, \"H\": 0.0}, \"family\": \"e\", \"familySubtype\": null, \"magneticCircuit\": "
            "\"open\", \"name\": \"E 32/16/11\", \"type\": \"standard\"}, \"type\": \"two-piece set\"}, "
            "\"manufacturerInfo\": null, \"name\": \"My Core\"}"));
        auto winding = OpenMagnetics::CoilWrapper(
            json::parse("{\"bobbin\": \"Dummy\", \"functionalDescription\": [{\"isolationSide\": \"primary\", \"name\": "
                        "\"Primary\", \"numberParallels\": 1, \"numberTurns\": 33, \"wire\": \"Dummy\"}], "
                        "\"layersDescription\": null, \"sectionsDescription\": null, \"turnsDescription\": null}"));
        auto operatingPoint = OpenMagnetics::OperatingPoint(
            json::parse("{\"conditions\": {\"ambientRelativeHumidity\": null, \"ambientTemperature\": 37.0, "
                        "\"cooling\": null, \"name\": null}, \"excitationsPerWinding\": [{\"current\": {\"harmonics\": "
                        "null, \"processed\": null, \"waveform\": {\"ancillaryLabel\": null, \"data\": [-5.0, 5.0, "
                        "-5.0], \"numberPeriods\": null, \"time\": [0.0, 2.4999999999999998e-06, 1e-05]}}, "
                        "\"frequency\": 100000.0, \"magneticFieldStrength\": null, \"magneticFluxDensity\": null, "
                        "\"magnetizingCurrent\": null, \"name\": \"My Operating Point\"}], \"name\": null}"));

        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::string{models["gapReluctance"]});

        OpenMagnetics::InputsWrapper::scale_time_to_frequency(operatingPoint, 215684);

        OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];

        auto magneticFluxDensity =
            magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint).second;

        excitation.set_magnetic_flux_density(magneticFluxDensity);
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();

        auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
        auto coreLosses = coreLossesModel->get_core_losses(core, excitation, temperature);


        auto frequencyFromCoreLosses = coreLossesModel->get_frequency_from_core_losses(core, magneticFluxDensity, temperature, coreLosses.get_core_losses()); 
        CHECK_CLOSE(excitation.get_frequency(), frequencyFromCoreLosses, frequencyFromCoreLosses * maxError);

    }
    TEST(Frequency_From_Losses_Micrometals) {
        auto models =
            json::parse(R"({"coreLosses": "PROPRIETARY", "coreTemperature": "MANIKTALA", "gapReluctance": "ZHANG"})");
        auto core = OpenMagnetics::CoreWrapper(json::parse(
            R"({"functionalDescription": {"gapping": [{"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 0.001, "sectionDimensions": null, "shape": null, "type": "subtractive"}, {"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}, {"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}], "material": "MS 40", "numberStacks": 1, "shape": {"aliases": [], "dimensions": {"A": 0.0391, "B": 0.0198, "C": 0.0125, "D": 0.0146, "E": 0.030100000000000002, "F": 0.0125, "G": 0.0, "H": 0.0}, "family": "etd", "familySubtype": "1", "magneticCircuit": null, "name": "ETD 39/20/13", "type": "standard"}, "type": "two-piece set"}, "manufacturerInfo": null, "name": "My Core"})"));
        auto winding = OpenMagnetics::CoilWrapper(json::parse(
            R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary", "numberParallels": 1, "numberTurns": 43, "wire": "Dummy"}], "layersDescription": null, "sectionsDescription": null, "turnsDescription": null})"));
        auto operatingPoint = OpenMagnetics::OperatingPoint(json::parse(
            R"({"conditions": {"ambientRelativeHumidity": null, "ambientTemperature": 25.0, "cooling": null, "name": null}, "excitationsPerWinding": [{"frequency": 100000.0, "magneticFieldStrength": null, "magneticFluxDensity": null, "magnetizingCurrent": null, "name": "My Operating Point", "voltage": {"harmonics": null, "processed": null, "waveform": {"ancillaryLabel": null, "data": [688.5, 688.5, -229.49999999999995, -229.49999999999995, 688.5], "numberPeriods": null, "time": [0.0, 2.4999999999999998e-06, 2.4999999999999998e-06, 1e-05, 1e-05]}}}], "name": null})"));
        OpenMagnetics::InputsWrapper::scale_time_to_frequency(operatingPoint, 123987);

        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::string{models["gapReluctance"]});

        OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];

        auto magneticFluxDensity =
            magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint).second;

        excitation.set_magnetic_flux_density(magneticFluxDensity);
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();

        auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
        auto coreLosses = coreLossesModel->get_core_losses(core, excitation, temperature);

        auto frequencyFromCoreLosses = coreLossesModel->get_frequency_from_core_losses(core, magneticFluxDensity, temperature, coreLosses.get_core_losses()); 
        CHECK_CLOSE(excitation.get_frequency(), frequencyFromCoreLosses, frequencyFromCoreLosses * maxError);

    }
}
SUITE(MagneticFluxDensityFromCoreLosses) {
    double maxError = 0.05;
    TEST(MagneticFluxDensity_From_Losses_Steinmetz) {
        auto models = json::parse("{\"coreLosses\": \"STEINMETZ\", \"gapReluctance\": \"BALAKRISHNAN\"}");
        auto core = OpenMagnetics::CoreWrapper(json::parse(
            "{\"functionalDescription\": {\"gapping\": [{\"area\": 9.8e-05, \"coordinates\": [0.0, 0.0001, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011301, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 0.0002, \"sectionDimensions\": [0.0092, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"subtractive\"}, {\"area\": 4.7e-05, \"coordinates\": [0.0138, 0.0, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011498, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.004401, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"residual\"}, {\"area\": 4.7e-05, \"coordinates\": [-0.0138, 0.0, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011498, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.004401, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"residual\"}], \"material\": \"3C97\", \"numberStacks\": 1, \"shape\": {\"aliases\": [], \"dimensions\": "
            "{\"A\": 0.032, \"B\": 0.0161, \"C\": 0.01065, \"D\": 0.0115, \"E\": 0.0232, \"F\": 0.0092, \"G\": 0.0, "
            "\"H\": 0.0}, \"family\": \"e\", \"familySubtype\": null, \"magneticCircuit\": \"open\", \"name\": \"E "
            "32/16/11\", \"type\": \"standard\"}, \"type\": \"two-piece set\"}, \"manufacturerInfo\": null, \"name\": "
            "\"My Core\"}"));
        double temperature = 42;
        double frequency = 423568;
        OpenMagnetics::OperatingPointExcitation operatingPointExcitation;
        OpenMagnetics::SignalDescriptor magneticFluxDensity;
        OpenMagnetics::Processed processed;
        operatingPointExcitation.set_frequency(frequency);
        processed.set_label(OpenMagnetics::WaveformLabel::SINUSOIDAL);
        processed.set_offset(0);
        processed.set_peak(0.3);
        processed.set_peak_to_peak(0.6);
        magneticFluxDensity.set_processed(processed);
        operatingPointExcitation.set_magnetic_flux_density(magneticFluxDensity);

        auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
        auto coreLosses = coreLossesModel->get_core_losses(core, operatingPointExcitation, temperature);

        auto magneticFluxDensityFromCoreLosses = coreLossesModel->get_magnetic_flux_density_from_core_losses(core, operatingPointExcitation.get_frequency(), temperature, coreLosses.get_core_losses()); 
        double magneticFluxDensityFromCoreLossesPeak = magneticFluxDensityFromCoreLosses.get_processed().value().get_peak().value();
        CHECK_CLOSE(magneticFluxDensity.get_processed().value().get_peak().value(), magneticFluxDensityFromCoreLossesPeak, magneticFluxDensityFromCoreLossesPeak * maxError);
    }
    TEST(MagneticFluxDensity_From_Losses_Igse) {
        auto models = json::parse("{\"coreLosses\": \"IGSE\", \"gapReluctance\": \"BALAKRISHNAN\"}");
        auto core = OpenMagnetics::CoreWrapper(json::parse(
            "{\"functionalDescription\": {\"gapping\": [{\"area\": 9.8e-05, \"coordinates\": [0.0, 0.0001, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011301, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 0.0002, \"sectionDimensions\": [0.0092, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"subtractive\"}, {\"area\": 4.7e-05, \"coordinates\": [0.0138, 0.0, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011498, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.004401, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"residual\"}, {\"area\": 4.7e-05, \"coordinates\": [-0.0138, 0.0, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011498, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.004401, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"residual\"}], \"material\": \"3C97\", \"numberStacks\": 1, \"shape\": {\"aliases\": [], \"dimensions\": "
            "{\"A\": 0.032, \"B\": 0.0161, \"C\": 0.01065, \"D\": 0.0115, \"E\": 0.0232, \"F\": 0.0092, \"G\": 0.0, "
            "\"H\": 0.0}, \"family\": \"e\", \"familySubtype\": null, \"magneticCircuit\": \"open\", \"name\": \"E "
            "32/16/11\", \"type\": \"standard\"}, \"type\": \"two-piece set\"}, \"manufacturerInfo\": null, \"name\": "
            "\"My Core\"}"));
        double temperature = 42;
        double frequency = 423568;
        OpenMagnetics::OperatingPointExcitation operatingPointExcitation;
        OpenMagnetics::SignalDescriptor magneticFluxDensity;
        OpenMagnetics::Processed processed;
        operatingPointExcitation.set_frequency(frequency);
        processed.set_label(OpenMagnetics::WaveformLabel::SINUSOIDAL);
        processed.set_offset(0);
        processed.set_peak(0.3);
        processed.set_peak_to_peak(0.6);
        magneticFluxDensity.set_processed(processed);
        operatingPointExcitation.set_magnetic_flux_density(magneticFluxDensity);

        auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
        auto coreLosses = coreLossesModel->get_core_losses(core, operatingPointExcitation, temperature);

        auto magneticFluxDensityFromCoreLosses = coreLossesModel->get_magnetic_flux_density_from_core_losses(core, operatingPointExcitation.get_frequency(), temperature, coreLosses.get_core_losses()); 
        double magneticFluxDensityFromCoreLossesPeak = magneticFluxDensityFromCoreLosses.get_processed().value().get_peak().value();
        CHECK_CLOSE(magneticFluxDensity.get_processed().value().get_peak().value(), magneticFluxDensityFromCoreLossesPeak, magneticFluxDensityFromCoreLossesPeak * maxError);
    }

    TEST(MagneticFluxDensity_From_Losses_Magnetics) {
        auto models = json::parse("{\"coreLosses\": \"PROPRIETARY\", \"gapReluctance\": \"BALAKRISHNAN\"}");
        auto core = OpenMagnetics::CoreWrapper(json::parse(
            "{\"functionalDescription\": {\"gapping\": [{\"area\": 9.8e-05, \"coordinates\": [0.0, 0.0001, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011301, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 0.0002, \"sectionDimensions\": [0.0092, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"subtractive\"}, {\"area\": 4.7e-05, \"coordinates\": [0.0138, 0.0, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011498, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.004401, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"residual\"}, {\"area\": 4.7e-05, \"coordinates\": [-0.0138, 0.0, 0.0], "
            "\"distanceClosestNormalSurface\": 0.011498, \"distanceClosestParallelSurface\": 0.006999999999999999, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.004401, 0.01065], \"shape\": \"rectangular\", \"type\": "
            "\"residual\"}], \"material\": \"XFlux 19\", \"numberStacks\": 1, \"shape\": {\"aliases\": [], "
            "\"dimensions\": {\"A\": 0.032, \"B\": 0.0161, \"C\": 0.01065, \"D\": 0.0115, \"E\": 0.0232, \"F\": "
            "0.0092, \"G\": 0.0, \"H\": 0.0}, \"family\": \"e\", \"familySubtype\": null, \"magneticCircuit\": "
            "\"open\", \"name\": \"E 32/16/11\", \"type\": \"standard\"}, \"type\": \"two-piece set\"}, "
            "\"manufacturerInfo\": null, \"name\": \"My Core\"}"));
        auto winding = OpenMagnetics::CoilWrapper(
            json::parse("{\"bobbin\": \"Dummy\", \"functionalDescription\": [{\"isolationSide\": \"primary\", \"name\": "
                        "\"Primary\", \"numberParallels\": 1, \"numberTurns\": 33, \"wire\": \"Dummy\"}], "
                        "\"layersDescription\": null, \"sectionsDescription\": null, \"turnsDescription\": null}"));
        auto operatingPoint = OpenMagnetics::OperatingPoint(
            json::parse("{\"conditions\": {\"ambientRelativeHumidity\": null, \"ambientTemperature\": 37.0, "
                        "\"cooling\": null, \"name\": null}, \"excitationsPerWinding\": [{\"current\": {\"harmonics\": "
                        "null, \"processed\": null, \"waveform\": {\"ancillaryLabel\": null, \"data\": [-5.0, 5.0, "
                        "-5.0], \"numberPeriods\": null, \"time\": [0.0, 2.4999999999999998e-06, 1e-05]}}, "
                        "\"frequency\": 100000.0, \"magneticFieldStrength\": null, \"magneticFluxDensity\": null, "
                        "\"magnetizingCurrent\": null, \"name\": \"My Operating Point\"}], \"name\": null}"));

        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::string{models["gapReluctance"]});
        OpenMagnetics::InputsWrapper::scale_time_to_frequency(operatingPoint, 215684);

        OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];

        auto magneticFluxDensity =
            magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint).second;

        excitation.set_magnetic_flux_density(magneticFluxDensity);
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();

        auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
        auto coreLosses = coreLossesModel->get_core_losses(core, excitation, temperature);

        auto magneticFluxDensityFromCoreLosses = coreLossesModel->get_magnetic_flux_density_from_core_losses(core, excitation.get_frequency(), temperature, coreLosses.get_core_losses()); 
        double magneticFluxDensityFromCoreLossesPeak = magneticFluxDensityFromCoreLosses.get_processed().value().get_peak().value();
        CHECK_CLOSE(magneticFluxDensity.get_processed().value().get_peak().value(), magneticFluxDensityFromCoreLossesPeak, magneticFluxDensityFromCoreLossesPeak * maxError);
    }

    TEST(MagneticFluxDensity_From_Losses_Micrometals) {
        auto models =
            json::parse(R"({"coreLosses": "PROPRIETARY", "coreTemperature": "MANIKTALA", "gapReluctance": "ZHANG"})");
        auto core = OpenMagnetics::CoreWrapper(json::parse(
            R"({"functionalDescription": {"gapping": [{"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 0.001, "sectionDimensions": null, "shape": null, "type": "subtractive"}, {"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}, {"area": null, "coordinates": null, "distanceClosestNormalSurface": null, "distanceClosestParallelSurface": null, "length": 1e-05, "sectionDimensions": null, "shape": null, "type": "residual"}], "material": "MS 40", "numberStacks": 1, "shape": {"aliases": [], "dimensions": {"A": 0.0391, "B": 0.0198, "C": 0.0125, "D": 0.0146, "E": 0.030100000000000002, "F": 0.0125, "G": 0.0, "H": 0.0}, "family": "etd", "familySubtype": "1", "magneticCircuit": null, "name": "ETD 39/20/13", "type": "standard"}, "type": "two-piece set"}, "manufacturerInfo": null, "name": "My Core"})"));
        auto winding = OpenMagnetics::CoilWrapper(json::parse(
            R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary", "numberParallels": 1, "numberTurns": 43, "wire": "Dummy"}], "layersDescription": null, "sectionsDescription": null, "turnsDescription": null})"));
        auto operatingPoint = OpenMagnetics::OperatingPoint(json::parse(
            R"({"conditions": {"ambientRelativeHumidity": null, "ambientTemperature": 25.0, "cooling": null, "name": null}, "excitationsPerWinding": [{"frequency": 100000.0, "magneticFieldStrength": null, "magneticFluxDensity": null, "magnetizingCurrent": null, "name": "My Operating Point", "voltage": {"harmonics": null, "processed": null, "waveform": {"ancillaryLabel": null, "data": [688.5, 688.5, -229.49999999999995, -229.49999999999995, 688.5], "numberPeriods": null, "time": [0.0, 2.4999999999999998e-06, 2.4999999999999998e-06, 1e-05, 1e-05]}}}], "name": null})"));

        OpenMagnetics::InputsWrapper::scale_time_to_frequency(operatingPoint, 123987);
        OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::string{models["gapReluctance"]});

        OpenMagnetics::OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];

        auto magneticFluxDensity =
            magnetizing_inductance.calculate_inductance_and_magnetic_flux_density(core, winding, &operatingPoint).second;

        excitation.set_magnetic_flux_density(magneticFluxDensity);
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();

        auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
        auto coreLosses = coreLossesModel->get_core_losses(core, excitation, temperature);

        auto magneticFluxDensityFromCoreLosses = coreLossesModel->get_magnetic_flux_density_from_core_losses(core, excitation.get_frequency(), temperature, coreLosses.get_core_losses()); 
        double magneticFluxDensityFromCoreLossesPeak = magneticFluxDensityFromCoreLosses.get_processed().value().get_peak().value();
        CHECK_CLOSE(magneticFluxDensity.get_processed().value().get_peak().value(), magneticFluxDensityFromCoreLossesPeak, magneticFluxDensityFromCoreLossesPeak * maxError);
    }
}