#pragma once
#include "MAS.hpp"
#include "Constants.h"
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

enum class ResistivityModels : int {
    CORE_MATERIAL,
    WIRE_MATERIAL
};

using ResistivityMaterial = std::variant<CoreMaterial, WireMaterial>;


class ResistivityModel {
  private:
  protected:
  public:
    virtual double get_resistivity(ResistivityMaterial materialData, double temperature) = 0;
    static std::shared_ptr<ResistivityModel> factory(ResistivityModels modelName);
};

class ResistivityCoreMaterialModel : public ResistivityModel {
  public:
    double get_resistivity(ResistivityMaterial materialData, double temperature);
};

class ResistivityWireMaterialModel : public ResistivityModel {
  public:
    double get_resistivity(ResistivityMaterial materialData, double temperature);
};

} // namespace OpenMagnetics