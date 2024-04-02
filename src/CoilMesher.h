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


enum class CoilMesherModels : int {
    WANG,
    CENTER
};


class CoilMesher {
  private:
  protected:
    double _windingLossesHarmonicAmplitudeThreshold;
    double _quickModeForManyHarmonicsThreshold = 1;
  public:
    std::vector<Field> generate_mesh_inducing_coil(MagneticWrapper magnetic, OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold = Defaults().windingLossesHarmonicAmplitudeThreshold, std::optional<std::vector<int8_t>> customCurrentDirectionPerWinding = std::nullopt);
    std::vector<Field> generate_mesh_induced_coil(MagneticWrapper magnetic, OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold = Defaults().windingLossesHarmonicAmplitudeThreshold);
    std::vector<size_t> get_common_harmonic_indexes(CoilWrapper coil, OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold);
    static std::pair<Field, double> generate_mesh_induced_grid(MagneticWrapper magnetic, double frequency, size_t numberPointsX, size_t numberPointsY);
};

class CoilMesherModel {
  private:
  public:
    std::string method_name = "Default";
    virtual std::vector<FieldPoint> generate_mesh_inducing_turn(Turn turn, WireWrapper wire, std::optional<size_t> turnIndex, std::optional<double> turnLength, CoreWrapper core) = 0;
    virtual std::vector<FieldPoint> generate_mesh_induced_turn(Turn turn, WireWrapper wire, std::optional<size_t> turnIndex = std::nullopt) = 0;
    static std::shared_ptr<CoilMesherModel> factory(CoilMesherModels modelName);

};

class CoilMesherCenterModel : public CoilMesherModel {
  public:
    std::vector<FieldPoint> generate_mesh_inducing_turn(Turn turn, [[maybe_unused]]WireWrapper wire, std::optional<size_t> turnIndex, std::optional<double> turnLength, CoreWrapper core);
    std::vector<FieldPoint> generate_mesh_induced_turn(Turn turn, [[maybe_unused]]WireWrapper wire, std::optional<size_t> turnIndex = std::nullopt);
};

// // Based on Improved Analytical Calculation of High Frequency Winding Losses in Planar Inductors by Xiaohui Wang
// // https://sci-hub.wf/10.1109/ECCE.2018.8558397
class CoilMesherWangModel : public CoilMesherModel {
  public:
    std::vector<FieldPoint> generate_mesh_induced_turn(Turn turn, [[maybe_unused]]WireWrapper wire, std::optional<size_t> turnIndex = std::nullopt);
    std::vector<FieldPoint> generate_mesh_inducing_turn(Turn turn, [[maybe_unused]]WireWrapper wire, std::optional<size_t> turnIndex, std::optional<double> turnLength, [[maybe_unused]]CoreWrapper core);
};



} // namespace OpenMagnetics