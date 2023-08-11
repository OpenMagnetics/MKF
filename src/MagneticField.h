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
    static SignalDescriptor get_magnetic_flux(SignalDescriptor magnetizingCurrent,
                                                      double reluctance,
                                                      double numberTurns,
                                                      double frequency);
    static SignalDescriptor get_magnetic_flux_density(SignalDescriptor magneticFlux,
                                                              double area,
                                                              double frequency);
    static SignalDescriptor get_magnetic_field_strength(SignalDescriptor magneticFluxDensity,
                                                                double initialPermeability,
                                                                double frequency);
};

} // namespace OpenMagnetics