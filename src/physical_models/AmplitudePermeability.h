#pragma once
#include "Constants.h"
#include "Definitions.h"
#include "support/Utils.h"

#include "constructive_models/Core.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <map>
#include <numbers>
#include <streambuf>
#include <vector>

using namespace MAS;

namespace OpenMagnetics {

class AmplitudePermeability {
    private:
    protected:
    public:
        static std::optional<double> get_amplitude_permeability(CoreMaterial coreMaterial, std::optional<double> magneticFluxDensityPeak, std::optional<double> magneticFieldStrengthPeak, double temperature=defaults.ambientTemperature);
        static std::optional<double> get_amplitude_permeability(std::string coreMaterialName, std::optional<double> magneticFluxDensityPeak, std::optional<double> magneticFieldStrengthPeak, double temperature=defaults.ambientTemperature);

};

class BHLoopModel {
  private:
  protected:

  public:
    BHLoopModel() = default;
    virtual ~BHLoopModel() = default;
    virtual std::pair<Curve2D, Curve2D> get_hysteresis_loop(std::string coreMaterialName,
                                                            double temperature,
                                                            std::optional<double> magneticFluxDensityPeak,
                                                            std::optional<double> magneticFieldStrengthPeak) = 0;

    virtual std::pair<Curve2D, Curve2D> get_hysteresis_loop(CoreMaterial coreMaterial,
                                                            double temperature,
                                                            std::optional<double> magneticFluxDensityPeak,
                                                            std::optional<double> magneticFieldStrengthPeak) = 0;
};

// Based on Ferrite Core Loss for Power Magnetic Components Design and
// A Practical, Accurate and Very General Core Loss Model for Nonsinusoidal Waveforms by Waseem Roshen
// https://sci-hub.st/10.1109/20.278656
// https://sci-hub.st/10.1109/TPEL.2006.886608
class BHLoopRoshenModel : public BHLoopModel {
  private:
    std::string _modelName = "Roshen";
  public:
    std::pair<Curve2D, Curve2D> get_hysteresis_loop(std::string coreMaterialName,
                                                    double temperature,
                                                    std::optional<double> magneticFluxDensityPeak,
                                                    std::optional<double> magneticFieldStrengthPeak);
    std::pair<Curve2D, Curve2D> get_hysteresis_loop(CoreMaterial coreMaterial,
                                                    double temperature,
                                                    std::optional<double> magneticFluxDensityPeak,
                                                    std::optional<double> magneticFieldStrengthPeak);
    std::map<std::string, double> get_major_loop_parameters(double saturationMagneticFieldStrength,
                                                            double saturationMagneticFluxDensity,
                                                            double coerciveForce,
                                                            double remanence);

};

} // namespace OpenMagnetics