#pragma once
#include "Constants.h"
#include "CoreLosses.h"
#include "Defaults.h"
#include "CoreWrapper.h"
#include <cmath>
#include <MAS.hpp>


namespace OpenMagnetics {

class CoreCrossReferencer {
    public: 
        enum class CoreCrossReferencerFilters : int {
            PERMEANCE,
            CORE_LOSSES, 
            WINDING_WINDOW_AREA, 
            DIMENSIONS
        };
    protected:
        std::map<std::string, std::string> _models;
        std::string _log;
        std::optional<std::string> _onlyManufacturer;
        std::map<CoreCrossReferencerFilters, double> _weights;

        void logEntry(std::string entry) {
            _log += entry + "\n";
        }

    public:
        std::map<CoreCrossReferencerFilters, std::map<std::string, bool>> _filterConfiguration{
                { CoreCrossReferencerFilters::PERMEANCE,            { {"invert", true}, {"log", true} } },
                { CoreCrossReferencerFilters::CORE_LOSSES,          { {"invert", true}, {"log", true} } },
                { CoreCrossReferencerFilters::WINDING_WINDOW_AREA,  { {"invert", true}, {"log", true} } },
                { CoreCrossReferencerFilters::DIMENSIONS,           { {"invert", true}, {"log", true} } }
            };
        std::map<CoreCrossReferencerFilters, std::map<std::string, double>> _scorings;
        CoreCrossReferencer(std::map<std::string, std::string> models) {
            auto defaults = OpenMagnetics::Defaults();
            _models = models;
            if (models.find("gapReluctance") == models.end()) {
                _models["gapReluctance"] = magic_enum::enum_name(defaults.reluctanceModelDefault);
            }
            if (models.find("coreLosses") == models.end()) {
                _models["gapReluctance"] = magic_enum::enum_name(defaults.coreLossesModelDefault);
            }
            if (models.find("coreTemperature") == models.end()) {
                _models["gapReluctance"] = magic_enum::enum_name(defaults.coreTemperatureModelDefault);
            }

            _weights[CoreCrossReferencerFilters::PERMEANCE] = 1;
            _weights[CoreCrossReferencerFilters::CORE_LOSSES] = 0.5;
            _weights[CoreCrossReferencerFilters::WINDING_WINDOW_AREA] = 0.5;
            _weights[CoreCrossReferencerFilters::DIMENSIONS] = 0.1;
        }
        CoreCrossReferencer() {
            auto defaults = OpenMagnetics::Defaults();
            _models["gapReluctance"] = magic_enum::enum_name(defaults.reluctanceModelDefault);
            _models["coreLosses"] = magic_enum::enum_name(defaults.coreLossesModelDefault);
            _models["coreTemperature"] = magic_enum::enum_name(defaults.coreTemperatureModelDefault);

            _weights[CoreCrossReferencerFilters::PERMEANCE] = 1;
            _weights[CoreCrossReferencerFilters::CORE_LOSSES] = 0.5;
            _weights[CoreCrossReferencerFilters::WINDING_WINDOW_AREA] = 0.5;
            _weights[CoreCrossReferencerFilters::DIMENSIONS] = 0.1;
        }
        std::string read_log() {
            return _log;
        }
        std::map<std::string, std::map<CoreCrossReferencerFilters, double>> get_scorings(){
            return get_scorings(false);
        }

        std::map<std::string, std::map<CoreCrossReferencerFilters, double>> get_scorings(bool weighted);

        std::vector<std::pair<CoreWrapper, double>> get_cross_referenced_core(CoreWrapper referenceCore, int64_t referenceNumberTurns, InputsWrapper inputs, size_t maximumNumberResults=10);
        std::vector<std::pair<CoreWrapper, double>> get_cross_referenced_core(CoreWrapper referenceCore, int64_t referenceNumberTurns, InputsWrapper inputs, std::map<CoreCrossReferencerFilters, double> weights, size_t maximumNumberResults=10);
        std::vector<std::pair<CoreWrapper, double>> apply_filters(std::vector<std::pair<CoreWrapper, double>>* cores, CoreWrapper referenceCore, int64_t referenceNumberTurns, InputsWrapper inputs, std::map<CoreCrossReferencerFilters, double> weights, size_t maximumNumberResults, double limit);

        void use_only_manufacturer(std::string onlyManufacturer) {
            _onlyManufacturer = onlyManufacturer;
        }

    class MagneticCoreFilter {
        public:
            std::map<CoreCrossReferencerFilters, std::map<std::string, double>>* _scorings;
            std::map<CoreCrossReferencerFilters, std::map<std::string, bool>>* _filterConfiguration;

            void add_scoring(std::string name, CoreCrossReferencer::CoreCrossReferencerFilters filter, double scoring) {
                if (scoring != -1) {
                    (*_scorings)[filter][name] = scoring;
                }
            }
            void set_scorings(std::map<CoreCrossReferencerFilters, std::map<std::string, double>>* scorings) {
                _scorings = scorings;
            }
            void set_filter_configuration(std::map<CoreCrossReferencerFilters, std::map<std::string, bool>>* filterConfiguration) {
                _filterConfiguration = filterConfiguration;
            }
            MagneticCoreFilter(){
            }
            std::vector<std::pair<CoreWrapper, double>> filter_core(std::vector<CoreWrapper> unfilteredCores, CoreWrapper referenceCore, InputsWrapper inputs, double weight=1);
    };
    
    class MagneticCoreFilterPermeance : public MagneticCoreFilter {
        public:
            std::vector<std::pair<CoreWrapper, double>> filter_core(std::vector<std::pair<CoreWrapper, double>>* unfilteredCores, CoreWrapper referenceCore, InputsWrapper inputs, std::map<std::string, std::string> models, double weight=1, double limit=0.25);
    };
    
    class MagneticCoreFilterWindingWindowArea : public MagneticCoreFilter {
        public:
            std::vector<std::pair<CoreWrapper, double>> filter_core(std::vector<std::pair<CoreWrapper, double>>* unfilteredCores, CoreWrapper referenceCore, double weight=1, double limit=0.25);
    };
    
    class MagneticCoreFilterDimensions : public MagneticCoreFilter {
        public:
            std::vector<std::pair<CoreWrapper, double>> filter_core(std::vector<std::pair<CoreWrapper, double>>* unfilteredCores, CoreWrapper referenceCore, double weight=1, double limit=0.25);
    };

    class MagneticCoreFilterCoreLosses : public MagneticCoreFilter {
        private:
            std::vector<CoreLossesModels> _coreLossesModelNames;
            std::vector<std::pair<CoreLossesModels, std::shared_ptr<CoreLossesModel>>> _coreLossesModels;

            std::vector<double> _magneticFluxDensities = {0.01, 0.025, 0.05, 0.1, 0.2};
            std::vector<double> _frequencies = {20000, 50000, 100000, 250000, 500000};
        public:
            MagneticCoreFilterCoreLosses() {
                _coreLossesModelNames = {Defaults().coreLossesModelDefault, CoreLossesModels::PROPRIETARY, CoreLossesModels::IGSE, CoreLossesModels::ROSHEN};
                for (auto modelName : _coreLossesModelNames) {
                    _coreLossesModels.push_back(std::pair<CoreLossesModels, std::shared_ptr<CoreLossesModel>>{modelName, CoreLossesModel::factory(modelName)});
                }
            }
            std::vector<std::pair<CoreWrapper, double>> filter_core(std::vector<std::pair<CoreWrapper, double>>* unfilteredCores, CoreWrapper referenceCore, int64_t referenceNumberTurns, InputsWrapper inputs, std::map<std::string, std::string> models, double weight=1, double limit=0.25);
            std::pair<double, double> calculate_average_core_losses_and_magnetic_flux_density(CoreWrapper core, int64_t numberTurns, InputsWrapper inputs, std::map<std::string, std::string> models);
    };
};


} // namespace OpenMagnetics