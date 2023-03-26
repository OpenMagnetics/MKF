#pragma once
#include "Constants.h"

#include <CoreWrapper.h>
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

class InitialPermeability {
  private:
  protected:
  public:
    double get_initial_permeability(CoreMaterialDataOrNameUnion material,
                                    double* temperature = nullptr,
                                    double* magneticFieldDcBias = nullptr,
                                    double* frequency = nullptr);
};

} // namespace OpenMagnetics