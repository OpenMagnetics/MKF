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
            AREA_PRODUCT, 
            ENERGY_STORED, 
            COST, 
            EFFICIENCY,
            DIMENSIONS,
            MINIMUM_IMPEDANCE
        };
    protected:
        std::map<std::string, std::string> _models;
        std::string _log;
        bool _uniqueCoreShapes = false;
        std::map<CoreAdviserFilters, double> _weights;
        MagneticSimulator _magneticSimulator;
        WindingOhmicLosses _windingOhmicLosses;

        void logEntry(std::string entry);

    public:

        std::map<CoreAdviserFilters, std::map<std::string, bool>> _filterConfiguration{
                { CoreAdviserFilters::AREA_PRODUCT,          { {"invert", true}, {"log", true} } },
                { CoreAdviserFilters::ENERGY_STORED,         { {"invert", true}, {"log", true} } },
                { CoreAdviserFilters::COST,                  { {"invert", true}, {"log", true} } },
                { CoreAdviserFilters::EFFICIENCY,            { {"invert", true}, {"log", true} } },
                { CoreAdviserFilters::DIMENSIONS,            { {"invert", true}, {"log", true} } },
                { CoreAdviserFilters::MINIMUM_IMPEDANCE,     { {"invert", true}, {"log", true} } },
            };
        std::map<CoreAdviserFilters, std::map<std::string, double>> _scorings;
        std::map<CoreAdviserFilters, std::map<std::string, bool>> _validScorings;
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
        std::string read_log() {
            return _log;
        }
        std::map<std::string, std::map<CoreAdviserFilters, double>> get_scorings(){
            return get_scorings(false);
        }

        std::map<std::string, std::map<CoreAdviserFilters, double>> get_scorings(bool weighted);

        void set_unique_core_shapes(bool value) {
            _uniqueCoreShapes = value;
        }

        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, size_t maximumNumberResults=1);
        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, std::vector<Core>* cores, size_t maximumNumberResults=1);
        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumNumberResults=1);
        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, std::map<CoreAdviserFilters, double> weights, std::vector<Core>* cores, size_t maximumNumberResults=1);
        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, std::vector<Core>* cores, size_t maximumNumberResults, size_t maximumNumberCores);
        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, std::map<CoreAdviserFilters, double> weights, std::vector<Core>* cores, size_t maximumNumberResults, size_t maximumNumberCores);
        std::vector<std::pair<Mas, double>> get_advised_core(Inputs inputs, std::vector<CoreShape>* shapes, std::vector<CoreMaterial>* materials, size_t maximumNumberResults=1);

        Mas post_process_core(Magnetic magnetic, Inputs inputs);
        std::vector<std::pair<Mas, double>> apply_filters(std::vector<std::pair<Magnetic, double>>* magnetics, Inputs inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults);
        std::vector<std::pair<Mas, double>> apply_fixed_filters(std::vector<std::pair<Magnetic, double>>* magnetics, Inputs inputs, size_t maximumNumberResults);
        std::vector<std::pair<Magnetic, double>> create_magnetic_dataset(Inputs inputs, std::vector<Core>* cores, bool includeStacks, bool onlyMaterialsForFilters=false);
        std::vector<std::pair<Magnetic, double>> create_magnetic_dataset(Inputs inputs, std::vector<CoreShape>* shapes, bool includeStacks, bool onlyMaterialsForFilters=false);
        void expand_magnetic_dataset_with_stacks(Inputs inputs, std::vector<Core>* cores, std::vector<std::pair<Magnetic, double>>* magnetics);
        void add_gapping(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs);
    
    class MagneticCoreFilter {
        public:
            std::map<CoreAdviserFilters, std::map<std::string, double>>* _scorings;
            std::map<CoreAdviserFilters, std::map<std::string, bool>>* _validScorings;
            std::map<CoreAdviserFilters, std::map<std::string, bool>>* _filterConfiguration;

            void add_scoring(std::string name, CoreAdviser::CoreAdviserFilters filter, double scoring, bool firstFilter) {
                if (firstFilter) {
                    (*_validScorings)[filter][name] = true;
                }
                if (scoring != -1) {
                    (*_scorings)[filter][name] = scoring;
                }
            }
            void set_scorings(std::map<CoreAdviserFilters, std::map<std::string, double>>* scorings) {
                _scorings = scorings;
            }
            void set_valid_scorings(std::map<CoreAdviserFilters, std::map<std::string, bool>>* validScorings) {
                _validScorings = validScorings;
            }
            void set_filter_configuration(std::map<CoreAdviserFilters, std::map<std::string, bool>>* filterConfiguration) {
                _filterConfiguration = filterConfiguration;
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

};


} // namespace OpenMagnetics