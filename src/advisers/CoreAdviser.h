#pragma once
#include "advisers/MagneticFilter.h"
#include "physical_models/CoreLosses.h"
#include "physical_models/Impedance.h"
#include "physical_models/MagneticEnergy.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "constructive_models/Coil.h"
#include "processors/MagneticSimulator.h"
#include "processors/Outputs.h"
#include "processors/Inputs.h"
#include "constructive_models/Core.h"
#include "constructive_models/Mas.h"
#include <cmath>
#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {

class CoreAdviser {
    public: 
        enum class CoreAdviserFilters : int {
            COST,
            EFFICIENCY,
            DIMENSIONS
        };

        enum class CoreAdviserModes : int {
            AVAILABLE_CORES,
            STANDARD_CORES,
            CUSTOM_CORES
        };
    protected:
        std::map<std::string, std::string> _models;
        bool _uniqueCoreShapes = false;
        std::map<CoreAdviserFilters, double> _weights;
        MagneticSimulator _magneticSimulator;
        WindingOhmicLosses _windingOhmicLosses;
        Application _application = Application::POWER;
        CoreAdviserModes _mode = CoreAdviserModes::AVAILABLE_CORES;


    public:

        std::map<CoreAdviserFilters, std::map<std::string, bool>> _filterConfiguration{
                { CoreAdviserFilters::COST,                  { {"invert", true}, {"log", true} } },
                { CoreAdviserFilters::EFFICIENCY,            { {"invert", true}, {"log", true} } },
                { CoreAdviserFilters::DIMENSIONS,            { {"invert", true}, {"log", true} } }
            };
        std::map<CoreAdviserFilters, std::map<std::string, double>> _scorings;
        CoreAdviser(std::map<std::string, std::string> models) {
            auto defaults = OpenMagnetics::Defaults();
            _models = models;
            if (models.find("gapReluctance") == models.end()) {
                _models["gapReluctance"] = magic_enum::enum_name(defaults.reluctanceModelDefault);
            }
        }
        CoreAdviser() {
            auto defaults = OpenMagnetics::Defaults();
            _models["gapReluctance"] = magic_enum::enum_name(defaults.reluctanceModelDefault);
            _models["coreLosses"] = magic_enum::enum_name(defaults.coreLossesModelDefault);
            _models["coreTemperature"] = magic_enum::enum_name(defaults.coreTemperatureModelDefault);
        }
        std::map<std::string, std::map<CoreAdviserFilters, double>> get_scorings(bool weighted = false);

        void set_unique_core_shapes(bool value);
        bool get_unique_core_shapes();
        void set_application(Application value);
        Application get_application();
        void set_mode(CoreAdviserModes value);
        CoreAdviserModes get_mode();

        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, size_t maximumNumberResults=1);
        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumNumberResults=1);
        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, std::vector<Core>* cores, size_t maximumNumberResults=1);
        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, std::map<CoreAdviserFilters, double> weights, std::vector<Core>* cores, size_t maximumNumberResults=1);
        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, std::vector<Core>* cores, size_t maximumNumberResults, size_t maximumNumberCores);
        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, std::map<CoreAdviserFilters, double> weights, std::vector<Core>* cores, size_t maximumNumberResults, size_t maximumNumberCores);
        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, std::vector<CoreShape>* shapes, size_t maximumNumberResults=1);

        Mas post_process_core(Magnetic magnetic, Inputs inputs);
        std::vector<std::pair<Mas, double>> filter_available_cores_power_application(std::vector<std::pair<Magnetic, double>>* magnetics, Inputs inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults);
        std::vector<std::pair<Mas, double>> filter_available_cores_suppression_application(std::vector<std::pair<Magnetic, double>>* magnetics, Inputs inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults);
        std::vector<std::pair<Mas, double>> filter_standard_cores_power_application(std::vector<std::pair<Magnetic, double>>* magnetics, Inputs inputs, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults);
        std::vector<std::pair<Mas, double>> filter_standard_cores_interference_suppression_application(std::vector<std::pair<Magnetic, double>>* magnetics, Inputs inputs, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults);
        std::vector<std::pair<Magnetic, double>> create_magnetic_dataset(Inputs inputs, std::vector<Core>* cores, bool includeStacks);
        std::vector<std::pair<Magnetic, double>> create_magnetic_dataset(Inputs inputs, std::vector<CoreShape>* shapes, bool includeStacks);
        void expand_magnetic_dataset_with_stacks(Inputs inputs, std::vector<Core>* cores, std::vector<std::pair<Magnetic, double>>* magnetics);
        std::vector<std::pair<Magnetic, double>> add_powder_materials(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs);
        std::vector<std::pair<Magnetic, double>> add_ferrite_materials_by_losses(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs);
        std::vector<std::pair<Magnetic, double>> add_ferrite_materials_by_impedance(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs);
        bool should_include_powder(Inputs inputs);
    
    class MagneticCoreFilter {
        public:
            std::map<CoreAdviserFilters, std::map<std::string, double>>* _scorings;
            std::map<CoreAdviserFilters, std::map<std::string, bool>>* _filterConfiguration;
            bool _useCache = true;

            void add_scoring(std::string name, CoreAdviser::CoreAdviserFilters filter, double scoring) {
                if (scoring != -1) {
                    if ((*_scorings)[filter].find(name) == (*_scorings)[filter].end()) {
                        (*_scorings)[filter][name] = 0;
                    }
                    auto oldScoring = (*_scorings)[filter][name];
                    (*_scorings)[filter][name] = oldScoring + scoring;
                }
            }
            void set_scorings(std::map<CoreAdviserFilters, std::map<std::string, double>>* scorings) {
                _scorings = scorings;
            }
            void set_filter_configuration(std::map<CoreAdviserFilters, std::map<std::string, bool>>* filterConfiguration) {
                _filterConfiguration = filterConfiguration;
            }
            void set_cache_usage(bool status) {
                _useCache = status;
            }
            MagneticCoreFilter(){
            }
            std::vector<std::pair<Magnetic, double>> filter_magnetics(std::vector<Magnetic> unfilteredMagnetics, Inputs inputs, double weight=1);
    };
    
    class MagneticCoreFilterAreaProduct : public MagneticCoreFilter {
        private:
            MagneticFilterAreaProduct _filter;

        public:
            MagneticCoreFilterAreaProduct(Inputs inputs);
            std::vector<std::pair<Magnetic, double>> filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight=1, bool firstFilter=false);
    };
    
    class MagneticCoreFilterEnergyStored : public MagneticCoreFilter {
        private:
            MagneticFilterEnergyStored _filter;

        public:
            MagneticCoreFilterEnergyStored(Inputs inputs, std::map<std::string, std::string> models);
            std::vector<std::pair<Magnetic, double>> filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight=1, bool firstFilter=false);
    };
    
    class MagneticCoreFilterFringingFactor : public MagneticCoreFilter {
        private:
            MagneticFilterFringingFactor _filter;

        public:
            MagneticCoreFilterFringingFactor(Inputs inputs, std::map<std::string, std::string> models);
            std::vector<std::pair<Magnetic, double>> filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight=1, bool firstFilter=false);
    };
    
    class MagneticCoreFilterCost : public MagneticCoreFilter {
        private:
            MagneticFilterEstimatedCost _filter;

        public:
            MagneticCoreFilterCost(Inputs inputs);
            std::vector<std::pair<Magnetic, double>> filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight=1, bool firstFilter=false);
    };
    
    class MagneticCoreFilterLosses : public MagneticCoreFilter {
        private:
            MagneticFilterCoreAndDcLosses _filter;

        public:
            MagneticCoreFilterLosses(Inputs inputs, std::map<std::string, std::string> models);
            std::vector<std::pair<Magnetic, double>> filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight=1, bool firstFilter=false);
    };
    
    class MagneticCoreFilterDimensions : public MagneticCoreFilter {
        private:
            MagneticFilterDimensions _filter;
        public:
            std::vector<std::pair<Magnetic, double>> filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, double weight=1, bool firstFilter=false);
    };
    
    class MagneticCoreFilterMinimumImpedance : public MagneticCoreFilter {
        private:
            MagneticFilterCoreMinimumImpedance _filter;
        public:
            std::vector<std::pair<Magnetic, double>> filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight=1, bool firstFilter=false);
    };
    
    class MagneticCoreFilterMagneticInductance : public MagneticCoreFilter {
        private:
            MagneticFilterMagnetizingInductance _filter;
        public:
            std::vector<std::pair<Magnetic, double>> filter_magnetics(std::vector<std::pair<Magnetic, double>>* unfilteredMagnetics, Inputs inputs, double weight=1, bool firstFilter=false);
    };

};

void from_json(const json & j, CoreAdviser::CoreAdviserModes & x);
void to_json(json & j, const CoreAdviser::CoreAdviserModes & x);

inline void from_json(const json & j, CoreAdviser::CoreAdviserModes & x) {
    if (j == "available cores") x = CoreAdviser::CoreAdviserModes::AVAILABLE_CORES;
    else if (j == "standard cores") x = CoreAdviser::CoreAdviserModes::STANDARD_CORES;
    else if (j == "custom cores") x = CoreAdviser::CoreAdviserModes::CUSTOM_CORES;
    else { throw std::runtime_error("Input JSON does not conform to schema!"); }
}

inline void to_json(json & j, const CoreAdviser::CoreAdviserModes & x) {
    switch (x) {
        case CoreAdviser::CoreAdviserModes::AVAILABLE_CORES: j = "available cores"; break;
        case CoreAdviser::CoreAdviserModes::STANDARD_CORES: j = "standard cores"; break;
        case CoreAdviser::CoreAdviserModes::CUSTOM_CORES: j = "custom cores"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"CoreAdviser::CoreAdviserModes\": " + std::to_string(static_cast<int>(x)));
    }
}

} // namespace OpenMagnetics