#pragma once
#include "Constants.h"

#include <InputsWrapper.h>
#include <MAS.hpp>
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

class MagneticField {
  private:
  protected:
  public:
    static ElectromagneticParameter get_magnetic_flux(ElectromagneticParameter magnetizingCurrent,
                                                      double reluctance,
                                                      double numberTurns,
                                                      double frequency);
    static ElectromagneticParameter get_magnetic_flux_density(ElectromagneticParameter magneticFlux,
                                                              double area,
                                                              double frequency);
    static ElectromagneticParameter get_magnetic_field_strength(ElectromagneticParameter magneticFluxDensity,
                                                                double initialPermeability,
                                                                double frequency);
};

} // namespace OpenMagnetics