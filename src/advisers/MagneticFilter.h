#pragma once
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
        std::shared_ptr<CoreLossesModel> _coreLossesModelSteinmetz = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "STEINMETZ"}}));
        std::shared_ptr<CoreLossesModel> _coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "PROPRIETARY"}}));
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
        std::shared_ptr<CoreLossesModel> _coreLossesModelSteinmetz = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "STEINMETZ"}}));
        std::shared_ptr<CoreLossesModel> _coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "PROPRIETARY"}}));
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
        std::pair<bool, double> evaluate_magnetic(CoilFunctionalDescription winding, Section section);
};

class MagneticFilterAreaWithParallels : public MagneticFilter {
    public:
        MagneticFilterAreaWithParallels() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        std::pair<bool, double> evaluate_magnetic(CoilFunctionalDescription winding, Section section, double numberSections, double sectionArea, bool allowNotFit);
};

class MagneticFilterEffectiveResistance : public MagneticFilter {
    private:
        double _maximumCurrent;
    public:
        MagneticFilterEffectiveResistance() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        std::pair<bool, double> evaluate_magnetic(CoilFunctionalDescription winding, double effectivefrequency, double temperature);
};

class MagneticFilterProximityFactor : public MagneticFilter {
    private:
        double _maximumCurrent;
    public:
        MagneticFilterProximityFactor() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        std::pair<bool, double> evaluate_magnetic(CoilFunctionalDescription winding, double effectiveSkinDepth, double temperature);
};

class MagneticFilterSolidInsulationRequirements : public MagneticFilter {
    private:
        double _maximumCurrent;
    public:
        MagneticFilterSolidInsulationRequirements() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        std::pair<bool, double> evaluate_magnetic(CoilFunctionalDescription winding, WireSolidInsulationRequirements wireSolidInsulationRequirements);
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

class MagneticFilterSkinLossesDensity : public MagneticFilter {
    public:
        MagneticFilterSkinLossesDensity() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        std::pair<bool, double> evaluate_magnetic(CoilFunctionalDescription winding, SignalDescriptor current, double temperature);
};

class MagneticFilterFringingFactor : public MagneticFilter {
    private:
        std::map<std::string, std::string> _models;
        MagneticEnergy _magneticEnergy;
        double _requiredMagneticEnergy;
        double _fringingFactorLitmit = 1.2;
        std::shared_ptr<ReluctanceModel> _reluctanceModel;

    public:
        MagneticFilterFringingFactor() {};
        MagneticFilterFringingFactor(Inputs inputs);
        MagneticFilterFringingFactor(Inputs inputs, std::map<std::string, std::string> models);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
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


// // Nice to have in the future
// class MagneticFilterMaximumWeight : public MagneticFilter {
//     public:
//         MagneticFilterMaximumWeight() {};
//         std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
// };

// // Nice to have in the future
// class MagneticFilterTemperature : public MagneticFilter {
//     public:
//         MagneticFilterTemperature() {};
//         std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
// };

// // Nice to have in the future
// class MagneticFilterStrayCapacitance : public MagneticFilter {
//     public:
//         MagneticFilterStrayCapacitance() {};
//         std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
// };

} // namespace OpenMagnetics