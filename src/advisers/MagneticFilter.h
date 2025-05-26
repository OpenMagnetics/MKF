#pragma once
#include "Constants.h"
#include "Defaults.h"
#include "physical_models/CoreLosses.h"
#include "physical_models/Impedance.h"
#include "physical_models/MagneticEnergy.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "constructive_models/Coil.h"
#include "processors/Outputs.h"
#include "processors/Inputs.h"
#include "constructive_models/Core.h"
#include "constructive_models/Mas.h"
#include <cmath>
#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {

class MagneticFilter {
    public: 
        enum class MagneticFilters : int {
            AREA_PRODUCT, 
            ENERGY_STORED, 
            COST, 
            EFFICIENCY,
            DIMENSIONS,
            MINIMUM_IMPEDANCE
        };

    public:

        std::map<MagneticFilters, std::map<std::string, bool>> _filterConfiguration{
                { MagneticFilters::AREA_PRODUCT,          { {"invert", true}, {"log", true} } },
                { MagneticFilters::ENERGY_STORED,         { {"invert", true}, {"log", true} } },
                { MagneticFilters::COST,                  { {"invert", true}, {"log", true} } },
                { MagneticFilters::EFFICIENCY,            { {"invert", true}, {"log", true} } },
                { MagneticFilters::DIMENSIONS,            { {"invert", true}, {"log", true} } },
                { MagneticFilters::MINIMUM_IMPEDANCE,     { {"invert", true}, {"log", true} } },
            };
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

    public:
        MagneticFilterEnergyStored() {};
        MagneticFilterEnergyStored(std::map<std::string, std::string> models);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
};

class MagneticFilterCost : public MagneticFilter {
    private:
        double _estimatedParallels = 0;
        double _estimatedWireTotalArea = 0;
        double _wireAirFillingFactor = 0;
        double _skinDepth = 0;
        double _averageMarginInWindingWindow = 0;

    public:
        MagneticFilterCost() {};
        MagneticFilterCost(Inputs inputs);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
};

class MagneticFilterLosses : public MagneticFilter {
    private:
        MagnetizingInductance _magnetizingInductance;
        WindingOhmicLosses _windingOhmicLosses;
        std::map<std::string, std::string> _models;
        std::shared_ptr<CoreLossesModel> _coreLossesModelSteinmetz;
        std::shared_ptr<CoreLossesModel> _coreLossesModelProprietary;
        double _maximumPowerMean = 0;

    public:
        MagneticFilterLosses() {};
        MagneticFilterLosses(Inputs inputs, std::map<std::string, std::string> models);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
};

class MagneticFilterDimensions : public MagneticFilter {
    public:
        MagneticFilterDimensions() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
};

class MagneticFilterMinimumImpedance : public MagneticFilter {
    private:
        Impedance _impedanceModel;
    public:
        MagneticFilterMinimumImpedance() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs);
};

class MagneticFilterAreaNoParallels : public MagneticFilter {
    private:
        int _maximumNumberParallels;
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

} // namespace OpenMagnetics