#pragma once
#include "physical_models/Temperature.h"
#include "physical_models/CoreLosses.h"
#include "physical_models/Impedance.h"
#include "physical_models/MagneticEnergy.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "constructive_models/Coil.h"
#include "processors/Inputs.h"
#include "processors/MagneticSimulator.h"
#include "processors/Outputs.h"
#include "constructive_models/Core.h"
#include "constructive_models/Mas.h"
#include "Definitions.h"
#include <cmath>
#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {

class MagneticFilter {
    public: 
        static std::shared_ptr<MagneticFilter> factory(MagneticFilters filterName, std::optional<Inputs> inputs = std::nullopt);

        MagneticFilter() { };
        virtual ~MagneticFilter() = default;
        virtual std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr) = 0;
};

class MagneticFilterAreaProduct : public MagneticFilter {
    private:
        std::map<std::string, double> _materialScaledMagneticFluxDensities;
        std::map<std::string, double> _bobbinFillingFactors;
        OperatingPointExcitation _operatingPointExcitation;
        std::vector<double> _areaProductRequiredPreCalculations;
        WindingSkinEffectLosses _windingSkinEffectLossesModel;
        std::shared_ptr<CoreLossesModel> _coreLossesModelSteinmetz;
        std::shared_ptr<CoreLossesModel> _coreLossesModelProprietary;
        double _averageMarginInWindingWindow = 0;
        double _magneticFluxDensityReference = 0.18;

    public:
        MagneticFilterAreaProduct() {};
        MagneticFilterAreaProduct(Inputs inputs);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        double get_estimated_area_product_required(Inputs inputs);
};

class MagneticFilterEnergyStored : public MagneticFilter {
    private:
        std::map<std::string, std::string> _models;
        MagneticEnergy _magneticEnergy;
        double _requiredMagneticEnergy;

    public:
        MagneticFilterEnergyStored() {};
        MagneticFilterEnergyStored(Inputs inputs);
        MagneticFilterEnergyStored(Inputs inputs, std::map<std::string, std::string> models);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterCost : public MagneticFilter {
    public:
        MagneticFilterCost() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterEstimatedCost : public MagneticFilter {
    private:
        double _estimatedParallels = 0;
        double _estimatedWireTotalArea = 0;
        double _wireAirFillingFactor = 0;
        double _skinDepth = 0;
        double _averageMarginInWindingWindow = 0;

    public:
        MagneticFilterEstimatedCost() {};
        MagneticFilterEstimatedCost(Inputs inputs);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterCoreAndDcLosses : public MagneticFilter {
    private:
        MagnetizingInductance _magnetizingInductance;
        WindingOhmicLosses _windingOhmicLosses;
        std::map<std::string, std::string> _models;
        std::shared_ptr<CoreLossesModel> _coreLossesModelSteinmetz = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Steinmetz"}}));
        std::shared_ptr<CoreLossesModel> _coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}));
        double _maximumPowerMean = 0;
    public:

        MagneticFilterCoreAndDcLosses();
        MagneticFilterCoreAndDcLosses(Inputs inputs);
        MagneticFilterCoreAndDcLosses(Inputs inputs, std::map<std::string, std::string> models);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterCoreDcAndSkinLosses : public MagneticFilter {
    private:
        MagnetizingInductance _magnetizingInductance;
        WindingOhmicLosses _windingOhmicLosses;
        WindingSkinEffectLosses _windingSkinEffectLosses;
        std::map<std::string, std::string> _models;
        std::shared_ptr<CoreLossesModel> _coreLossesModelSteinmetz = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Steinmetz"}}));
        std::shared_ptr<CoreLossesModel> _coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}));
        double _maximumPowerMean = 0;
    public:

        MagneticFilterCoreDcAndSkinLosses();
        MagneticFilterCoreDcAndSkinLosses(Inputs inputs);
        MagneticFilterCoreDcAndSkinLosses(Inputs inputs, std::map<std::string, std::string> models);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterLosses : public MagneticFilter {
    private:
        std::map<std::string, std::string> _models;
        MagneticSimulator _magneticSimulator;
    public:
        MagneticFilterLosses() {};
        MagneticFilterLosses(std::map<std::string, std::string> models);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterLossesNoProximity : public MagneticFilter {
    private:
        std::map<std::string, std::string> _models;
        WindingOhmicLosses _windingOhmicLosses;
        WindingSkinEffectLosses _windingSkinEffectLosses;
        MagneticSimulator _magneticSimulator;
    public:
        MagneticFilterLossesNoProximity() {};
        MagneticFilterLossesNoProximity(std::map<std::string, std::string> models);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterDimensions : public MagneticFilter {
    public:
        MagneticFilterDimensions() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

/**
 * @class MagneticFilterTurnCount
 * @brief Manufacturability / copper-burden proxy: turns × core size.
 *
 * Score = (Σ N_i) × core_width   (lower is better, caller inverts).
 *
 * The core-size factor is mandatory. For a fixed |Z| target on a CMC,
 * Z = µ_complex · ω · N² · Aₑ / lₑ, so N ∝ 1/√(µ·Aₑ/lₑ): bigger cores
 * minimise N. Scoring N alone therefore monotonically rewards the
 * largest core in the catalogue — which is the opposite of "manufacturable"
 * (more wire per turn, more copper, more cost). Multiplying by
 * core_width (OD for toroids, A-dim for two-piece sets) gives a simple
 * 1-D wire-length proxy that lets a small T 25 with N=20 outrank a big
 * T 140 with N=13. See TestTopologyCmc::Test_Cmc_AdviserMustNotPickOversizedToroid_WizardDefaults.
 *
 * Read-only: requires N to be already populated on the coil's
 * functional_description (e.g. via add_initial_turns_by_inductance /
 * filterMinimumImpedance). Returns valid=true with score=0 if no turns
 * have been assigned yet (caller-side gate ensures ordering).
 */
class MagneticFilterTurnCount : public MagneticFilter {
    public:
        MagneticFilterTurnCount() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterCoreMinimumImpedance : public MagneticFilter {
    private:
        Impedance _impedanceModel;
    public:
        MagneticFilterCoreMinimumImpedance() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterAreaNoParallels : public MagneticFilter {
    private:
        int _maximumNumberParallels = defaults.maximumNumberParallels;
    public:
        MagneticFilterAreaNoParallels() {};
        MagneticFilterAreaNoParallels(int maximumNumberParallels);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        std::pair<bool, double> evaluate_magnetic(Winding winding, Section section);
};

class MagneticFilterAreaWithParallels : public MagneticFilter {
    public:
        MagneticFilterAreaWithParallels() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        std::pair<bool, double> evaluate_magnetic(Winding winding, Section section, double numberSections, double sectionArea, bool allowNotFit);
};

class MagneticFilterEffectiveResistance : public MagneticFilter {
    private:
        double _maximumCurrent;
    public:
        MagneticFilterEffectiveResistance() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        std::pair<bool, double> evaluate_magnetic(Winding winding, double effectivefrequency, double temperature);
};

class MagneticFilterProximityFactor : public MagneticFilter {
    private:
        double _maximumCurrent;
    public:
        MagneticFilterProximityFactor() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        std::pair<bool, double> evaluate_magnetic(Winding winding, double effectiveSkinDepth, double temperature);
};

class MagneticFilterSolidInsulationRequirements : public MagneticFilter {
    private:
        double _maximumCurrent;
    public:
        MagneticFilterSolidInsulationRequirements() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        std::pair<bool, double> evaluate_magnetic(Winding winding, WireSolidInsulationRequirements wireSolidInsulationRequirements);
};

class MagneticFilterTurnsRatios : public MagneticFilter {
    public:
        MagneticFilterTurnsRatios() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterMaximumDimensions : public MagneticFilter {
    public:
        MagneticFilterMaximumDimensions() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterSaturation : public MagneticFilter {
    public:
        MagneticFilterSaturation() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterDcCurrentDensity : public MagneticFilter {
    public:
        MagneticFilterDcCurrentDensity() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterEffectiveCurrentDensity : public MagneticFilter {
    public:
        MagneticFilterEffectiveCurrentDensity() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterImpedance : public MagneticFilter {
    public:
        MagneticFilterImpedance() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterMagnetizingInductance : public MagneticFilter {
    public:
        MagneticFilterMagnetizingInductance() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

/**
 * @class MagneticFilterLeakageInductance
 * @brief Filter for evaluating and scoring magnetic designs based on leakage inductance.
 *
 * For Common Mode Chokes (CMCs), lower leakage inductance indicates tighter coupling
 * between windings, which is essential for effective common-mode rejection.
 * The coupling coefficient k = 1 - (Lk/Lm) should be close to 1.
 *
 * For transformers requiring low leakage (e.g., LLC resonant converters),
 * this filter validates against the leakage_inductance design requirement.
 *
 * @note Returns leakage inductance in Henries. Lower values score better when inverted.
 */
class MagneticFilterLeakageInductance : public MagneticFilter {
    public:
        MagneticFilterLeakageInductance() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterSkinLossesDensity : public MagneticFilter {
    public:
        MagneticFilterSkinLossesDensity() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        std::pair<bool, double> evaluate_magnetic(Winding winding, SignalDescriptor current, double temperature);
};

class MagneticFilterFringingFactor : public MagneticFilter {
    private:
        std::map<std::string, std::string> _models;
        MagneticEnergy _magneticEnergy;
        double _requiredMagneticEnergy;
        double _fringingFactorLitmit = 1.2;
        std::shared_ptr<ReluctanceModel> _reluctanceModel;
        // OPTIMIZATION: Cache for fringing factor calculations to avoid redundant computation
        static std::map<std::string, double> _fringingFactorCache;

    public:
        MagneticFilterFringingFactor() {};
        MagneticFilterFringingFactor(Inputs inputs);
        MagneticFilterFringingFactor(Inputs inputs, std::map<std::string, std::string> models);
        void set_fringing_factor_limit(double limit) { _fringingFactorLitmit = limit; }
        double get_fringing_factor_limit() const { return _fringingFactorLitmit; }
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        // OPTIMIZATION: Clear cache when needed (e.g., between different design runs)
        static void clear_cache();
};

class MagneticFilterVolume : public MagneticFilter {
    public:
        MagneticFilterVolume() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterArea : public MagneticFilter {
    public:
        MagneticFilterArea() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterHeight : public MagneticFilter {
    public:
        MagneticFilterHeight() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterTemperatureRise : public MagneticFilter {
    private:
        MagneticFilterLossesNoProximity _magneticFilterLossesNoProximity;
    public:
        MagneticFilterTemperatureRise() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterLossesTimesVolume : public MagneticFilter {
    public:
        MagneticFilterLossesTimesVolume() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterVolumeTimesTemperatureRise : public MagneticFilter {
    private:
        MagneticFilterTemperatureRise _magneticFilterTemperatureRise;
    public:
        MagneticFilterVolumeTimesTemperatureRise() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterLossesTimesVolumeTimesTemperatureRise : public MagneticFilter {
    private:
        MagneticFilterTemperatureRise _magneticFilterTemperatureRise;
    public:
        MagneticFilterLossesTimesVolumeTimesTemperatureRise() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterLossesNoProximityTimesVolume : public MagneticFilter {
    private:
        MagneticFilterLossesNoProximity _magneticFilterLossesNoProximity;
    public:
        MagneticFilterLossesNoProximityTimesVolume() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterLossesNoProximityTimesVolumeTimesTemperatureRise : public MagneticFilter {
    private:
        MagneticFilterTemperatureRise _magneticFilterTemperatureRise;
        MagneticFilterLossesNoProximity _magneticFilterLossesNoProximity;
    public:
        MagneticFilterLossesNoProximityTimesVolumeTimesTemperatureRise() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagnetomotiveForce : public MagneticFilter {
    public:
        MagnetomotiveForce() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};


// // Nice to have in the future
// class MagneticFilterMaximumWeight : public MagneticFilter {
//     public:
//         MagneticFilterMaximumWeight() {};
//         std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
// };

class MagneticFilterTemperature : public MagneticFilter {
    double _maximumTemperature = 130.0;
    // Orchestrator, not a fixed model: selects per material from its available
    // volumetric-losses methods (Steinmetz family, proprietary, loss factor)
    CoreLosses _coreLosses;
    MagnetizingInductance _magnetizingInductance;
public:
    MagneticFilterTemperature() {};
    MagneticFilterTemperature(Inputs inputs, double maximumTemperature);
    std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs,
                                              std::vector<Outputs>* outputs = nullptr);
};

// // Nice to have in the future
// class MagneticFilterStrayCapacitance : public MagneticFilter {
//     public:
//         MagneticFilterStrayCapacitance() {};
//         std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
// };

} // namespace OpenMagnetics