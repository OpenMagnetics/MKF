#pragma once
#include "Constants.h"
#include "MAS.hpp"
#include "Resistivity.h"
#include "CoilWrapper.h"
#include "WireWrapper.h"
#include "Utils.h"
#include "Models.h"
#include "Defaults.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <map>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {


class WindingSkinEffectLossesModel {
    private:
    protected:
        std::map<size_t, std::map<double, std::map<double, double>>> _skinFactorPerWirePerFrequencyPerTemperature;
        std::optional<double> try_get_skin_factor(WireWrapper wire,  double frequency, double temperature);
        void set_skin_factor(WireWrapper wire, double frequency, double temperature, double skinFactor);

    public:
        std::string methodName = "Default";
        virtual double calculate_turn_losses(WireWrapper wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0) = 0;
        static std::shared_ptr<WindingSkinEffectLossesModel> factory(WindingSkinEffectLossesModels modelName);
};


class WindingSkinEffectLosses {
    private:
    protected:
    public:
      static std::shared_ptr<WindingSkinEffectLossesModel> get_model(WireType wireType);
      static double calculate_skin_depth(WireMaterialDataOrNameUnion material, double frequency, double temperature);
      static double calculate_skin_depth(WireWrapper wire, double frequency, double temperature);
      static WindingLossesOutput calculate_skin_effect_losses(CoilWrapper coil, double temperature, WindingLossesOutput windingLossesOutput, double windingLossesHarmonicAmplitudeThreshold = Defaults().windingLossesHarmonicAmplitudeThreshold);
      static std::pair<double, std::vector<std::pair<double, double>>> calculate_skin_effect_losses_per_meter(WireWrapper wire, SignalDescriptor current, double temperature, double currentDivider = 1, double windingLossesHarmonicAmplitudeThreshold = Defaults().windingLossesHarmonicAmplitudeThreshold);

};

// // Based on Effects of eddy currents in transformer windings by P.L. Dowell
// // http://www.cpdee.ufmg.br/~troliveira/docs/aulas/fontes/05247417.pdf
// class WindingSkinEffectLossesDowellModel : public WindingSkinEffectLossesModel {
//   public:
//     static double get_skin_effect_losses_per_meter(Turn turn);
//     static std::map<std::string, double> get_skin_effect_losses_per_winding(Winding winding);
// };

// Based on Winding Resistance and Power Loss of Inductors With Litz and Solid-Round Wires by Rafal P. Wojda
// https://sci-hub.wf/https://ieeexplore.ieee.org/document/8329131
class WindingSkinEffectLossesWojdaModel : public WindingSkinEffectLossesModel {
  public:
    std::string methodName = "Wojda";
    double calculate_skin_factor(const WireWrapper& wire, double frequency, double temperature);
    double calculate_penetration_ratio(const WireWrapper& wire, double frequency, double temperature);
    double calculate_turn_losses(WireWrapper wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0);

};

// Based on Induktivitäten in der Leistungselektronik: Spulen, Trafos und ihre parasitären Eigenschaften by Manfred Albach
// https://libgen.rocks/get.php?md5=94b7f2906f53602f19892d7f1dabd929&key=YMKCEJOWB653PYLL
// https://sci-hub.wf/10.1109/tpel.2011.2143729
class WindingSkinEffectLossesAlbachModel : public WindingSkinEffectLossesModel {
  public:
    double calculate_skin_factor(const WireWrapper& wire, double frequency, double temperature);
    double calculate_turn_losses(WireWrapper wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0);
};

// Based on The Ac Resistance Of Rectangular Conductors by Alan Payne
// https://www.researchgate.net/publication/351307928_THE_AC_RESISTANCE_OF_RECTANGULAR_CONDUCTORS
class WindingSkinEffectLossesPayneModel : public WindingSkinEffectLossesModel {
  public:
    double calculate_turn_losses(WireWrapper wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0);
};

// // Based on An Improved Calculation of Proximity-Effect Loss in High-Frequency Windings of Round Conductors by Xi Nan
// // http://inductor.thayerschool.org/papers/newcalc.pdf
// class WindingSkinEffectLossesNanModel : public WindingSkinEffectLossesModel {
//   public:
//     static double get_skin_effect_losses_per_meter(Turn turn);
//     static std::map<std::string, double> get_skin_effect_losses_per_winding(Winding winding);
// };

// // Based on An Analytical Correction to Dowell's Equation for Inductor and Transformer Winding Losses Using Cylindrical Coordinates by Marian K. Kazimierczuk
// // https://sci-hub.st/https://ieeexplore.ieee.org/document/8665954
// class WindingSkinEffectLossesKazimierczukModel : public WindingSkinEffectLossesModel {
//   public:
//     static double get_skin_effect_losses_per_meter(Turn turn);
//     static std::map<std::string, double> get_skin_effect_losses_per_winding(Winding winding);
// };

// Based on A Simple Technique to Evaluate Winding Losses Including Two-Dimensional Edge Effects by Nasser H. Kutkut
// https://sci-hub.wf/10.1109/63.712319
class WindingSkinEffectLossesKutkutModel : public WindingSkinEffectLossesModel {
  public:
    double calculate_turn_losses(WireWrapper wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0);
};

// Based on Two-Dimensional Skin Effect in Power Foils for High-Frequency Applications by Ashraf W. Lotfi
// https://sci-hub.wf/https://ieeexplore.ieee.org/document/364775
class WindingSkinEffectLossesLotfiModel : public WindingSkinEffectLossesModel {
  public:
    double calculate_turn_losses(WireWrapper wire, [[maybe_unused]] double dcLossTurn, double frequency, double temperature, double currentRms = 0);
};

// Based on A New Approach to Analyse Conduction Losses in High Frequency Magnetic Components by J.A. Ferreira
// https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9485268
class WindingSkinEffectLossesFerreiraModel : public WindingSkinEffectLossesModel {
  public:
    double calculate_skin_factor(const WireWrapper& wire, double frequency, double temperature);
    double calculate_turn_losses(WireWrapper wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0);
};

// // Based on A New Model for the Determination of Copper Losses in Transformer Windings with Arbitrary Conductor Distribution under High Frequency Sinusoidal Excitation by G. S. Dimitrakakis
// // https://sci-hub.wf/10.1109/epe.2007.4417574
// class WindingSkinEffectLossesDimitrakakisModel : public WindingSkinEffectLossesModel {
//   public:
//     static double get_skin_effect_losses_per_meter(Turn turn);
//     static std::map<std::string, double> get_skin_effect_losses_per_winding(Winding winding);
// };

// // Based on Modeling And Multi-objective Optimization Of Inductive Power Components by Jonas Mühlethaler
// // https://www.pes-publications.ee.ethz.ch/uploads/tx_ethpublications/Diss_Muehlethaler.pdf
// class WindingSkinEffectLossesMuehlethalerModel : public WindingSkinEffectLossesModel {
//   public:
//     static double get_skin_effect_losses_per_meter(Turn turn);
//     static std::map<std::string, double> get_skin_effect_losses_per_winding(Winding winding);
// };

// // Based on Improved Analytical Calculation of High Frequency Winding Losses in Planar Inductors by Xiaohui Wang
// // https://sci-hub.wf/10.1109/ECCE.2018.8558397
// class WindingSkinEffectLossesWangModel : public WindingSkinEffectLossesModel {
//   public:
//     static double get_skin_effect_losses_per_meter(Turn turn);
//     static std::map<std::string, double> get_skin_effect_losses_per_winding(Winding winding);
// };

// // Based on Power Losses Calculations in Windings of Gapped Magnetic Components by Fermín A. Holguin
// // https://oa.upm.es/37148/1/INVE_MEM_2014_197924.pdf
// class WindingSkinEffectLossesHolguinModel : public WindingSkinEffectLossesModel {
//   public:
//     static double get_skin_effect_losses_per_meter(Turn turn);
//     static std::map<std::string, double> get_skin_effect_losses_per_winding(Winding winding);
// };


// // Based on Multiple Layer Series Connected Winding Design For Minimum Losses by M.P. Perry
// // https://sci-hub.wf/10.1109/tpas.1979.319520
// class WindingSkinEffectLossesPerryModel : public WindingSkinEffectLossesModel {
//   public:
//     static double get_skin_effect_losses_per_meter(Turn turn);
//     static std::map<std::string, double> get_skin_effect_losses_per_winding(Winding winding);
// };


} // namespace OpenMagnetics