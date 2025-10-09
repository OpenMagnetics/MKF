#pragma once
#include "Constants.h"
#include "physical_models/CoreLosses.h"
#include "Defaults.h"
#include "constructive_models/Core.h"
#include <cmath>
#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {

class CoreCrossReferencer {
    public: 
    protected:
        std::map<std::string, std::string> _models;
        std::string _log;
        std::optional<std::string> _onlyManufacturer;
        bool _onlyReferenceMaterial = false;
        double _limit = 1;
        std::map<CoreCrossReferencerFilters, double> _weights;

    public:
        std::map<CoreCrossReferencerFilters, std::map<std::string, bool>> _filterConfiguration{
                { CoreCrossReferencerFilters::PERMEANCE,            { {"invert", true}, {"log", false} } },
                { CoreCrossReferencerFilters::CORE_LOSSES,          { {"invert", true}, {"log", false} } },
                { CoreCrossReferencerFilters::SATURATION,           { {"invert", true}, {"log", false} } },
                { CoreCrossReferencerFilters::WINDING_WINDOW_AREA,  { {"invert", true}, {"log", false} } },
                { CoreCrossReferencerFilters::EFFECTIVE_AREA,       { {"invert", true}, {"log", false} } },
                { CoreCrossReferencerFilters::ENVELOPING_VOLUME,    { {"invert", true}, {"log", false} } }
            };
        std::map<CoreCrossReferencerFilters, std::map<std::string, double>> _scorings;
        std::map<CoreCrossReferencerFilters, std::map<std::string, bool>> _validScorings;
        std::map<CoreCrossReferencerFilters, std::map<std::string, double>> _scoredValues;
        CoreCrossReferencer(std::map<std::string, std::string> models) {
            auto defaults = OpenMagnetics::Defaults();
            _models = models;
            if (models.find("gapReluctance") == models.end()) {
                _models["gapReluctance"] = to_string(defaults.reluctanceModelDefault);
            }
            if (models.find("coreLosses") == models.end()) {
                _models["coreLosses"] = to_string(defaults.coreLossesModelDefault);
            }
            if (models.find("coreTemperature") == models.end()) {
                _models["coreTemperature"] = to_string(defaults.coreTemperatureModelDefault);
            }

            _weights[CoreCrossReferencerFilters::PERMEANCE] = 1;
            _weights[CoreCrossReferencerFilters::SATURATION] = 0.5;
            _weights[CoreCrossReferencerFilters::CORE_LOSSES] = 0.5;
            _weights[CoreCrossReferencerFilters::EFFECTIVE_AREA] = 0.5;
            _weights[CoreCrossReferencerFilters::WINDING_WINDOW_AREA] = 0.5;
            _weights[CoreCrossReferencerFilters::ENVELOPING_VOLUME] = 0.1;
        }
        CoreCrossReferencer() {
            auto defaults = OpenMagnetics::Defaults();
            _models["gapReluctance"] = to_string(defaults.reluctanceModelDefault);
            _models["coreLosses"] = to_string(defaults.coreLossesModelDefault);
            _models["coreTemperature"] = to_string(defaults.coreTemperatureModelDefault);

            _weights[CoreCrossReferencerFilters::PERMEANCE] = 1;
            _weights[CoreCrossReferencerFilters::SATURATION] = 0.5;
            _weights[CoreCrossReferencerFilters::CORE_LOSSES] = 0.5;
            _weights[CoreCrossReferencerFilters::EFFECTIVE_AREA] = 0.5;
            _weights[CoreCrossReferencerFilters::WINDING_WINDOW_AREA] = 0.5;
            _weights[CoreCrossReferencerFilters::ENVELOPING_VOLUME] = 0.1;
        }
        std::string read_log() {
            return _log;
        }
        std::map<std::string, std::map<CoreCrossReferencerFilters, double>> get_scorings(){
            return get_scorings(false);
        }

        std::map<std::string, std::map<CoreCrossReferencerFilters, double>> get_scored_values();
        std::map<std::string, std::map<CoreCrossReferencerFilters, double>> get_scorings(bool weighted);

        std::vector<std::pair<Core, double>> get_cross_referenced_core(Core referenceCore, int64_t referenceNumberTurns, Inputs inputs, size_t maximumNumberResults=10);
        std::vector<std::pair<Core, double>> get_cross_referenced_core(Core referenceCore, int64_t referenceNumberTurns, Inputs inputs, std::map<CoreCrossReferencerFilters, double> weights, size_t maximumNumberResults=10);
        std::vector<std::pair<Core, double>> apply_filters(std::vector<std::pair<Core, double>>* cores, Core referenceCore, int64_t referenceNumberTurns, Inputs inputs, std::map<CoreCrossReferencerFilters, double> weights, size_t maximumNumberResults, double limit);

        void use_only_manufacturer(std::string onlyManufacturer) {
            _onlyManufacturer = onlyManufacturer;
        }
        void use_only_reference_material(bool value) {
            _onlyReferenceMaterial = value;
        }
        void set_limit(double value) {
            _limit = value;
        }

    class MagneticCoreFilter {
        public:
            std::map<CoreCrossReferencerFilters, std::map<std::string, double>>* _scorings;
            std::map<CoreCrossReferencerFilters, std::map<std::string, bool>>* _validScorings;
            std::map<CoreCrossReferencerFilters, std::map<std::string, double>>* _scoredValues;
            std::map<CoreCrossReferencerFilters, std::map<std::string, bool>>* _filterConfiguration;

            void add_scoring(std::string name, CoreCrossReferencerFilters filter, double scoring) {
                if (std::isnan(scoring)) {
                    throw std::invalid_argument("scoring cannot be nan");
                }

                if (scoring != -1) {
                    (*_validScorings)[filter][name] = true;
                    (*_scorings)[filter][name] = scoring;
                }
            }
            void add_scored_value(std::string name, CoreCrossReferencerFilters filter, double scoredValues) {
                if (scoredValues != -1) {
                    (*_scoredValues)[filter][name] = scoredValues;
                }
            }
            void set_scorings(std::map<CoreCrossReferencerFilters, std::map<std::string, double>>* scorings) {
                _scorings = scorings;
            }
            void set_valid_scorings(std::map<CoreCrossReferencerFilters, std::map<std::string, bool>>* validScorings) {
                _validScorings = validScorings;
            }
            void set_scored_value(std::map<CoreCrossReferencerFilters, std::map<std::string, double>>* scoredValues) {
                _scoredValues = scoredValues;
            }
            void set_filter_configuration(std::map<CoreCrossReferencerFilters, std::map<std::string, bool>>* filterConfiguration) {
                _filterConfiguration = filterConfiguration;
            }
            MagneticCoreFilter(){
            }
            std::vector<std::pair<Core, double>> filter_core(std::vector<Core> unfilteredCores, Core referenceCore, Inputs inputs, double weight=1);
    };
    
    class MagneticCoreFilterPermeance : public MagneticCoreFilter {
        public:
            std::vector<std::pair<Core, double>> filter_core(std::vector<std::pair<Core, double>>* unfilteredCores, Core referenceCore, Inputs inputs, std::map<std::string, std::string> models, double weight=1, double limit=0.25);
    };
    
    class MagneticCoreFilterWindingWindowArea : public MagneticCoreFilter {
        public:
            std::vector<std::pair<Core, double>> filter_core(std::vector<std::pair<Core, double>>* unfilteredCores, Core referenceCore, double weight=1, double limit=0.25);
    };
    
    class MagneticCoreFilterEffectiveArea : public MagneticCoreFilter {
        public:
            std::vector<std::pair<Core, double>> filter_core(std::vector<std::pair<Core, double>>* unfilteredCores, Core referenceCore, double weight=1, double limit=0.25);
    };
    
    class MagneticCoreFilterEnvelopingVolume : public MagneticCoreFilter {
        public:
            std::vector<std::pair<Core, double>> filter_core(std::vector<std::pair<Core, double>>* unfilteredCores, Core referenceCore, double weight=1, double limit=0.25);
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
            std::vector<std::pair<Core, double>> filter_core(std::vector<std::pair<Core, double>>* unfilteredCores, Core referenceCore, int64_t referenceNumberTurns, Inputs inputs, std::map<std::string, std::string> models, double weight=1, double limit=0.25);
            std::pair<double, double> calculate_average_core_losses_and_magnetic_flux_density(Core core, int64_t numberTurns, Inputs inputs, std::map<std::string, std::string> models);
    };
};

} // namespace OpenMagnetics