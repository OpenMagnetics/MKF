#pragma once
#include "Constants.h"
#include "MAS.hpp"
#include "Resistivity.h"
#include "WindingWrapper.h"
#include "WireWrapper.h"
#include "Utils.h"

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

// enum class WindingSkinEffectLossesModels : int {
//     DOWELL,
//     WODJA,
//     ALBACH,
//     PAYNE,
//     // NAN,
//     VANDELAC_ZIOGAS,
//     KAZIMIERCZUK,
//     KUTKUT,
//     FERREIRA,
//     DIMITRAKAKIS,
//     WANG,
//     HOLGUIN,
//     PERRY,
// }

class WindingSkinEffectLossesModel {
  private:
  protected:
  public:
    // virtual double get_skin_effect_losses_per_meter(Turn turn) = 0;
    // virtual std::map<std::string, double> get_skin_effect_losses_per_winding(Winding winding) = 0;
    static double get_skin_depth(WireMaterialDataOrNameUnion material, double frequency, double temperature);
    // static double get_effective_current_density(Harmonics harmonics, WireWrapper wire, double temperature);

};

// // Based on Effects of eddy currents in transformer windings by P.L. Dowell
// // http://www.cpdee.ufmg.br/~troliveira/docs/aulas/fontes/05247417.pdf
// class WindingSkinEffectLossesDowellModel : public WindingSkinEffectLossesModel {
//   public:
//     static double get_skin_effect_losses_per_meter(Turn turn);
//     static std::map<std::string, double> get_skin_effect_losses_per_winding(Winding winding);
// };

// // Based on Winding Resistance and Power Loss of Inductors With Litz and Solid-Round Wires by Rafal P. Wojda
// // https://sci-hub.wf/https://ieeexplore.ieee.org/document/8329131
// class WindingSkinEffectLossesWodjaModel : public WindingSkinEffectLossesModel {
//   public:
//     static double get_skin_effect_losses_per_meter(Turn turn);
//     static std::map<std::string, double> get_skin_effect_losses_per_winding(Winding winding);
// };

// // Based on Induktivitäten in der Leistungselektronik: Spulen, Trafos und ihre parasitären Eigenschaften by Manfred Albach
// // https://libgen.rocks/get.php?md5=94b7f2906f53602f19892d7f1dabd929&key=YMKCEJOWB653PYLL
// class WindingSkinEffectLossesAlbachModel : public WindingSkinEffectLossesModel {
//   public:
//     static double get_skin_effect_losses_per_meter(Turn turn);
//     static std::map<std::string, double> get_skin_effect_losses_per_winding(Winding winding);
// };

// // Based on The Ac Resistance Of Rectangular Conductors by Alan Payne
// // https://www.researchgate.net/publication/351307928_THE_AC_RESISTANCE_OF_RECTANGULAR_CONDUCTORS
// class WindingSkinEffectLossesPayneModel : public WindingSkinEffectLossesModel {
//   public:
//     static double get_skin_effect_losses_per_meter(Turn turn);
//     static std::map<std::string, double> get_skin_effect_losses_per_winding(Winding winding);
// };

// // Based on An Improved Calculation of Proximity-Effect Loss in High-Frequency Windings of Round Conductors by Xi Nan
// // http://inductor.thayerschool.org/papers/newcalc.pdf
// class WindingSkinEffectLossesNanModel : public WindingSkinEffectLossesModel {
//   public:
//     static double get_skin_effect_losses_per_meter(Turn turn);
//     static std::map<std::string, double> get_skin_effect_losses_per_winding(Winding winding);
// };

// // Based on A Novel Approach for Minimizing High-Frequency Transformer Copper Losses by Jean-Pierre Vandelac And Phoivos D. Ziogas
// // https://sci-hub.wf/10.1109/63.17944
// class WindingSkinEffectLossesVandelacZiogasModel : public WindingSkinEffectLossesModel {
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

// // Based on A Simple Technique to Evaluate Winding Losses Including Two-Dimensional Edge Effects by Nasser H. Kutkut
// // https://sci-hub.wf/10.1109/63.712319
// class WindingSkinEffectLossesKutkutModel : public WindingSkinEffectLossesModel {
//   public:
//     static double get_skin_effect_losses_per_meter(Turn turn);
//     static std::map<std::string, double> get_skin_effect_losses_per_winding(Winding winding);
// };

// // Based on A New Approach to Analyse Conduction Losses in High Frequency Magnetic Components by J.A. Ferreira
// // https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9485268
// class WindingSkinEffectLossesFerreiraModel : public WindingSkinEffectLossesModel {
//   public:
//     static double get_skin_effect_losses_per_meter(Turn turn);
//     static std::map<std::string, double> get_skin_effect_losses_per_winding(Winding winding);
// };

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