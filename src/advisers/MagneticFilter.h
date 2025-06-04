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
        virtual std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs) = 0;
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

    public:
        MagneticFilterAreaProduct() {};
        MagneticFilterAreaProduct(Inputs inputs);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
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
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
};

class MagneticFilterCost : public MagneticFilter {
    public:
        MagneticFilterCost() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
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
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
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
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
};

class MagneticFilterLosses : public MagneticFilter {
    private:
        std::map<std::string, std::string> _models;
        MagneticSimulator _magneticSimulator;
    public:
        MagneticFilterLosses() {};
        MagneticFilterLosses(std::map<std::string, std::string> models);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
};

class MagneticFilterDimensions : public MagneticFilter {
    public:
        MagneticFilterDimensions() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
};

class MagneticFilterCoreMinimumImpedance : public MagneticFilter {
    private:
        Impedance _impedanceModel;
    public:
        MagneticFilterCoreMinimumImpedance() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
};

class MagneticFilterAreaNoParallels : public MagneticFilter {
    private:
        int _maximumNumberParallels = defaults.maximumNumberParallels;
    public:
        MagneticFilterAreaNoParallels() {};
        MagneticFilterAreaNoParallels(int maximumNumberParallels);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
        std::pair<bool, double> evaluate_magnetic(CoilFunctionalDescription winding, Section section);
};

class MagneticFilterAreaWithParallels : public MagneticFilter {
    public:
        MagneticFilterAreaWithParallels() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
        std::pair<bool, double> evaluate_magnetic(CoilFunctionalDescription winding, Section section, double numberSections, double sectionArea, bool allowNotFit);
};

class MagneticFilterEffectiveResistance : public MagneticFilter {
    private:
        double _maximumCurrent;
    public:
        MagneticFilterEffectiveResistance() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
        std::pair<bool, double> evaluate_magnetic(CoilFunctionalDescription winding, double effectivefrequency, double temperature);
};

class MagneticFilterProximityFactor : public MagneticFilter {
    private:
        double _maximumCurrent;
    public:
        MagneticFilterProximityFactor() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
        std::pair<bool, double> evaluate_magnetic(CoilFunctionalDescription winding, double effectiveSkinDepth, double temperature);
};

class MagneticFilterSolidInsulationRequirements : public MagneticFilter {
    private:
        double _maximumCurrent;
    public:
        MagneticFilterSolidInsulationRequirements() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
        std::pair<bool, double> evaluate_magnetic(CoilFunctionalDescription winding, WireSolidInsulationRequirements wireSolidInsulationRequirements);
};

class MagneticFilterTurnsRatios : public MagneticFilter {
    public:
        MagneticFilterTurnsRatios() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
};

class MagneticFilterMaximumDimensions : public MagneticFilter {
    public:
        MagneticFilterMaximumDimensions() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
};

class MagneticFilterSaturation : public MagneticFilter {
    public:
        MagneticFilterSaturation() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
};

class MagneticFilterDcCurrentDensity : public MagneticFilter {
    public:
        MagneticFilterDcCurrentDensity() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
};

class MagneticFilterEffectiveCurrentDensity : public MagneticFilter {
    public:
        MagneticFilterEffectiveCurrentDensity() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
};

class MagneticFilterImpedance : public MagneticFilter {
    public:
        MagneticFilterImpedance() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
};

class MagneticFilterMagnetizingInductance : public MagneticFilter {
    public:
        MagneticFilterMagnetizingInductance() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
};

// // Nice to have in the future
// class MagneticFilterMaximumWeight : public MagneticFilter {
//     public:
//         MagneticFilterMaximumWeight() {};
//         std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
// };

// // Nice to have in the future
// class MagneticFilterTemperature : public MagneticFilter {
//     public:
//         MagneticFilterTemperature() {};
//         std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
// };

// // Nice to have in the future
// class MagneticFilterStrayCapacitance : public MagneticFilter {
//     public:
//         MagneticFilterStrayCapacitance() {};
//         std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
// };

} // namespace OpenMagnetics