#pragma once
#include "Constants.h"
#include "MAS.hpp"
#include "CoilWrapper.h"
#include "WireWrapper.h"
#include "Utils.h"
#include "Models.h"

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


class WindingProximityEffectLossesModel {
  private:
  protected:
    std::map<size_t, std::map<double, std::map<double, double>>> _proximityFactorPerWirePerFrequencyPerTemperature;
  public:
    std::string methodName = "Default";
    virtual double calculate_turn_losses(WireWrapper wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature) = 0;
    std::optional<double> try_get_proximity_factor(WireWrapper wire, double frequency, double temperature);
    void set_proximity_factor(WireWrapper wire, double frequency, double temperature, double proximityFactor);
    static std::shared_ptr<WindingProximityEffectLossesModel> factory(WindingProximityEffectLossesModels modelName);
};

class WindingProximityEffectLosses {
  private:
  protected:
  public:
    static std::shared_ptr<WindingProximityEffectLossesModel> get_model(WireType wireType);
    static WindingLossesOutput calculate_proximity_effect_losses(CoilWrapper coil, double temperature, WindingLossesOutput windingLossesOutput, WindingWindowMagneticStrengthFieldOutput windingWindowMagneticStrengthFieldOutput);
    static std::pair<double, std::vector<std::pair<double, double>>> calculate_proximity_effect_losses_per_meter(WireWrapper wire, double temperature, std::vector<ComplexField> fields);

};

// Based on Measurement and Characterization of High Frequency Losses in Nonideal Litz Wires by Hans Rossmanith
// https://sci-hub.wf/10.1109/tpel.2011.2143729
class WindingProximityEffectLossesRossmanithModel : public WindingProximityEffectLossesModel {
  public:
    std::string methodName = "Rossmanith";
    double calculate_proximity_factor(WireWrapper wire, double frequency, double temperature);
    double calculate_turn_losses(WireWrapper wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature);

};

// Based on Improved Analytical Calculation of High Frequency Winding Losses in Planar Inductors by Xiaohui Wang
// https://sci-hub.wf/10.1109/ECCE.2018.8558397
class WindingProximityEffectLossesWangModel : public WindingProximityEffectLossesModel {
  public:
    std::string methodName = "Wang";
    double calculate_turn_losses(WireWrapper wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature);
};

// Based on A New Approach to Analyse Conduction Losses in High Frequency Magnetic Components by J.A. Ferreira
// https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9485268
class WindingProximityEffectLossesFerreiraModel : public WindingProximityEffectLossesModel {
  public:
    std::string methodName = "Ferreira";
    double calculate_proximity_factor(WireWrapper wire, double frequency, double temperature);
    double calculate_turn_losses(WireWrapper wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature);
};


// Based on Induktivitäten in der Leistungselektronik: Spulen, Trafos und ihre parasitären Eigenschaften by Manfred Albach
// https://libgen.rocks/get.php?md5=94b7f2906f53602f19892d7f1dabd929&key=YMKCEJOWB653PYLL
class WindingProximityEffectLossesAlbachModel : public WindingProximityEffectLossesModel {
  public:
    std::string methodName = "Albach";
    double calculate_turn_losses(WireWrapper wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature);
};

// Based on Eddy currents by Jiří Lammeraner
// https://archive.org/details/eddycurrents0000lamm
class WindingProximityEffectLossesLammeranerModel : public WindingProximityEffectLossesModel {
  public:
    std::string methodName = "Lammeraner";
    double calculate_proximity_factor(WireWrapper wire, double frequency, double temperature);
    double calculate_turn_losses(WireWrapper wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature);
};

} // namespace OpenMagnetics