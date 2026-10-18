#pragma once
#include "Constants.h"
#include "Definitions.h"
#include "support/Utils.h"

#include "constructive_models/Core.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <streambuf>
#include <vector>

using namespace MAS;

namespace OpenMagnetics {

enum class CoreLoopPropertyTier {
    MATERIAL,      // the material's own MAS value
    CLASS,         // the sourced members of its class (material type + composition)
    MATERIAL_TYPE  // the sourced members of its material type (proprietary or missing composition, or an empty class)
};

/**
 * @brief How a remanence or coercive force that MAS does not carry was estimated (ABT #1951).
 */
struct CoreLoopPropertyEstimate {
    CoreLoopPropertyTier tier;
    // The pool the estimate was drawn from, e.g. "class FERRITE MN_ZN",
    // "material type NANOCRYSTALLINE, round loop (initial permeability from 15000)".
    std::string pool;
    // Remanence: the pool's Br/Bs (the loop squareness of the class), applied to the material's own Bs.
    // Coercive force: the pool's Hc in A/m at the default ambient temperature.
    double classValue;
    std::vector<std::string> sourceMaterials;
    size_t numberDocuments = 0;
    size_t numberDataPoints = 0;
    std::string method;
};

/**
 * @brief A core material's remanence (T) or coercive force (A/m) at one temperature, and whether it
 * is the material's own value or an estimate.
 */
struct CoreLoopProperty {
    double value;
    CoreLoopPropertyTier tier;
    // Engaged if and only if MAS carries no value of this property for the material.
    std::optional<CoreLoopPropertyEstimate> estimate;
    bool is_estimated() const { return estimate.has_value(); }
};

/**
 * @brief THE resolvers of the remanence and coercive force the Roshen hysteresis loop is fitted
 * through (Alf, 2026-10-18, ABT #1951, approved like the heat-conductivity estimator of ABT #1752/#1753).
 *
 * The material's own MAS data always wins (Core::get_remanence / Core::get_coercive_force). Without
 * it, a flagged estimate from the materials of its class that carry a sourced value:
 * remanence = (median Br/Bs of the class) x the material's own Bs(T); coercive force = the median
 * Hc of the class at the default ambient temperature, going linearly to zero at the material's
 * Curie temperature past it (as tabulated data does, ABT #1736). Tiers, first that applies:
 * CLASS (material type + a physical composition), MATERIAL_TYPE. A NANOCRYSTALLINE material pools
 * only the members of its own loop family, flat loop (initial permeability below 15000) or round loop
 * (from 15000), never the other family's (Alf, 2026-10-19; derivation in AmplitudePermeability.cpp).
 * A value is sourced when a
 * non-retracted, non-derived provenance entry names the field and cites a document, or, for a value
 * no provenance entry names, when its record cites the maker's datasheet. Every number is taken
 * from MAS at call time; nothing is hardcoded. An estimate logs one WARNING per material and property.
 * Throws MaterialException(MATERIAL_DATA_MISSING) when neither the class nor the material type has a
 * sourced member (the class's loop is unknown), and MATERIAL_INVALID_PROPERTY for a pooled value
 * that is not a finite positive number.
 */
CoreLoopProperty resolve_core_material_remanence(const CoreMaterial& coreMaterial, double temperature);
CoreLoopProperty resolve_core_material_coercive_force(const CoreMaterial& coreMaterial, double temperature);

class AmplitudePermeability {
    private:
    protected:
    public:
        static std::optional<double> get_amplitude_permeability(const CoreMaterial& coreMaterial, std::optional<double> magneticFluxDensityPeak, std::optional<double> magneticFieldStrengthPeak, double temperature=defaults.ambientTemperature);
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

    virtual std::pair<Curve2D, Curve2D> get_hysteresis_loop(const CoreMaterial& coreMaterial,
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
    std::pair<Curve2D, Curve2D> get_hysteresis_loop(const CoreMaterial& coreMaterial,
                                                    double temperature,
                                                    std::optional<double> magneticFluxDensityPeak,
                                                    std::optional<double> magneticFieldStrengthPeak);
    std::map<std::string, double> get_major_loop_parameters(double saturationMagneticFieldStrength,
                                                            double saturationMagneticFluxDensity,
                                                            double coerciveForce,
                                                            double remanence);

};

} // namespace OpenMagnetics