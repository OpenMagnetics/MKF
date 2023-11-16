#pragma once
#include "Defaults.h"
#include "MAS.hpp"
#include "MagneticWrapper.h"
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


enum class CoilBreakerModels : int {
    SQUARE,
    RECTANGULAR,
    CENTER,
    WANG
};


class CoilBreaker {
  private:
  protected:
    double _windingLossesHarmonicAmplitudeThreshold;
    int _mirroringDimension = 0;
  public:
    WindingWindowCurrentFieldOutput breakdown_coil(MagneticWrapper magnetic, OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold = Defaults().windingLossesHarmonicAmplitudeThreshold);

    void set_mirroring_dimension(int mirroringDimension) {
        _mirroringDimension = mirroringDimension;
    }
};

class CoilBreakerModel {
  private:
  protected:
    int _mirroringDimension = 0;
  public:
    std::string method_name = "Default";
    virtual std::vector<FieldPoint> breakdown_turn(Turn turn, WireWrapper wire, double currentPeak, std::optional<size_t> turnIndex, std::optional<double> turnLength, CoreWrapper core) = 0;
    static std::shared_ptr<CoilBreakerModel> factory(CoilBreakerModels modelName);

    void set_mirroring_dimension(int mirroringDimension) {
        _mirroringDimension = mirroringDimension;
    }
};

class CoilBreakerCenterModel : public CoilBreakerModel {
  public:
    std::vector<FieldPoint> breakdown_turn(Turn turn, WireWrapper wire, double currentPeak, std::optional<size_t> turnIndex, std::optional<double> turnLength, CoreWrapper core);
};

// // Based on Improved Analytical Calculation of High Frequency Winding Losses in Planar Inductors by Xiaohui Wang
// // https://sci-hub.wf/10.1109/ECCE.2018.8558397
// class CoilBreakerWangModel : public CoilBreakerModel {
//   public:
//     static double get_skin_effect_losses_per_meter(Turn turn);
//     static std::map<std::string, double> get_skin_effect_losses_per_winding(Winding winding);
// };



} // namespace OpenMagnetics