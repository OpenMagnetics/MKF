#pragma once
#include "Constants.h"
#include "physical_models/CoreLosses.h"
#include "Defaults.h"
#include "constructive_models/Core.h"
#include <cmath>
#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {

class CoreMaterialCrossReferencer {
    public: 
        enum class CoreMaterialCrossReferencerFilters : int {
            INITIAL_PERMEABILITY,
            REMANENCE, 
            COERCIVE_FORCE, 
            SATURATION,
            CURIE_TEMPERATURE,
            VOLUMETRIC_LOSSES,
            RESISTIVITY
        };
    protected:
        std::map<std::string, std::string> _models;
        std::string _log;
        std::optional<std::string> _onlyManufacturer;
        std::map<CoreMaterialCrossReferencerFilters, double> _weights;


    public:
        std::map<CoreMaterialCrossReferencerFilters, std::map<std::string, bool>> _filterConfiguration{
                { CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY, { {"invert", true}, {"log", false} } },
                { CoreMaterialCrossReferencerFilters::REMANENCE,            { {"invert", true}, {"log", false} } },
                { CoreMaterialCrossReferencerFilters::COERCIVE_FORCE,       { {"invert", true}, {"log", false} } },
                { CoreMaterialCrossReferencerFilters::SATURATION,           { {"invert", true}, {"log", false} } },
                { CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE,    { {"invert", true}, {"log", false} } },
                { CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES,    { {"invert", true}, {"log", false} } },
                { CoreMaterialCrossReferencerFilters::RESISTIVITY,          { {"invert", true}, {"log", false} } },
            };
        std::map<CoreMaterialCrossReferencerFilters, std::map<std::string, double>> _scorings;
        std::map<CoreMaterialCrossReferencerFilters, std::map<std::string, double>> _scoredValues;
        CoreMaterialCrossReferencer(std::map<std::string, std::string> models) {
            auto defaults = OpenMagnetics::Defaults();
            _models = models;

            if (models.find("gapReluctance") == models.end()) {
                _models["gapReluctance"] = magic_enum::enum_name(defaults.reluctanceModelDefault);
            }
            if (models.find("coreLosses") == models.end()) {
                _models["coreLosses"] = magic_enum::enum_name(defaults.coreLossesModelDefault);
            }
            if (models.find("coreTemperature") == models.end()) {
                _models["coreTemperature"] = magic_enum::enum_name(defaults.coreTemperatureModelDefault);
            }

            _weights[CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY] = 0.5;
            _weights[CoreMaterialCrossReferencerFilters::REMANENCE] = 0.01;
            _weights[CoreMaterialCrossReferencerFilters::COERCIVE_FORCE] = 0.01;
            _weights[CoreMaterialCrossReferencerFilters::SATURATION] = 1;
            _weights[CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE] = 0.01;
            _weights[CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES] = 1;
            _weights[CoreMaterialCrossReferencerFilters::RESISTIVITY] = 0.2;
        }
        CoreMaterialCrossReferencer() {
            auto defaults = OpenMagnetics::Defaults();
            _models["gapReluctance"] = magic_enum::enum_name(defaults.reluctanceModelDefault);
            _models["coreLosses"] = magic_enum::enum_name(defaults.coreLossesModelDefault);
            _models["coreTemperature"] = magic_enum::enum_name(defaults.coreTemperatureModelDefault);

            _weights[CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY] = 0.5;
            _weights[CoreMaterialCrossReferencerFilters::REMANENCE] = 0.01;
            _weights[CoreMaterialCrossReferencerFilters::COERCIVE_FORCE] = 0.01;
            _weights[CoreMaterialCrossReferencerFilters::SATURATION] = 1;
            _weights[CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE] = 0.01;
            _weights[CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES] = 1;
            _weights[CoreMaterialCrossReferencerFilters::RESISTIVITY] = 0.2;
        }
        std::string read_log() {
            return _log;
        }
        std::map<std::string, std::map<CoreMaterialCrossReferencerFilters, double>> get_scorings(){
            return get_scorings(false);
        }

        std::map<std::string, std::map<CoreMaterialCrossReferencerFilters, double>> get_scored_values();
        std::map<std::string, std::map<CoreMaterialCrossReferencerFilters, double>> get_scorings(bool weighted);

        std::vector<std::pair<CoreMaterial, double>> get_cross_referenced_core_material(CoreMaterial referenceCoreMaterial, double temperature, size_t maximumNumberResults=10);
        std::vector<std::pair<CoreMaterial, double>> get_cross_referenced_core_material(CoreMaterial referenceCoreMaterial, double temperature, std::map<CoreMaterialCrossReferencerFilters, double> weights, size_t maximumNumberResults=10);
        std::vector<std::pair<CoreMaterial, double>> apply_filters(std::vector<std::pair<CoreMaterial, double>>* coreMaterials, CoreMaterial referenceCoreMaterial, double temperature, std::map<CoreMaterialCrossReferencerFilters, double> weights, size_t maximumNumberResults);

        void use_only_manufacturer(std::string onlyManufacturer) {
            _onlyManufacturer = onlyManufacturer;
        }

    class MagneticCoreFilter {
        public:
            std::map<CoreMaterialCrossReferencerFilters, std::map<std::string, double>>* _scorings;
            std::map<CoreMaterialCrossReferencerFilters, std::map<std::string, double>>* _scoredValues;
            std::map<CoreMaterialCrossReferencerFilters, std::map<std::string, bool>>* _filterConfiguration;

            void add_scoring(std::string name, CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters filter, double scoring) {
                if (std::isnan(scoring)) {
                    throw std::invalid_argument("scoring cannot be nan");
                }

                if (scoring != -1) {
                    (*_scorings)[filter][name] = scoring;
                }
            }
            void add_scored_value(std::string name, CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters filter, double scoredValues) {
                if (scoredValues != -1) {
                    (*_scoredValues)[filter][name] = scoredValues;
                }
            }
            void set_scorings(std::map<CoreMaterialCrossReferencerFilters, std::map<std::string, double>>* scorings) {
                _scorings = scorings;
            }
            void set_scored_value(std::map<CoreMaterialCrossReferencerFilters, std::map<std::string, double>>* scoredValues) {
                _scoredValues = scoredValues;
            }
            void set_filter_configuration(std::map<CoreMaterialCrossReferencerFilters, std::map<std::string, bool>>* filterConfiguration) {
                _filterConfiguration = filterConfiguration;
            }
            MagneticCoreFilter(){
            }
            std::vector<std::pair<CoreMaterial, double>> filter_core_materials(std::vector<CoreMaterial> unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, double weight=1);
    };
    
    class MagneticCoreFilterInitialPermeability : public MagneticCoreFilter {
        public:
            std::vector<std::pair<CoreMaterial, double>> filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, double weight=1);
    };
    
    class MagneticCoreFilterRemanence : public MagneticCoreFilter {
        public:
            std::vector<std::pair<CoreMaterial, double>> filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, double weight=1);
    };
    
    class MagneticCoreFilterCoerciveForce : public MagneticCoreFilter {
        public:
            std::vector<std::pair<CoreMaterial, double>> filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, double weight=1);
    };
    
    class MagneticCoreFilterSaturation : public MagneticCoreFilter {
        public:
            std::vector<std::pair<CoreMaterial, double>> filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, double weight=1);
    };
    
    class MagneticCoreFilterCurieTemperature : public MagneticCoreFilter {
        public:
            std::vector<std::pair<CoreMaterial, double>> filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, double weight=1);
    };
    
    class MagneticCoreFilterVolumetricLosses : public MagneticCoreFilter {
        private:
            std::vector<CoreLossesModels> _coreLossesModelNames;
            std::vector<std::pair<CoreLossesModels, std::shared_ptr<CoreLossesModel>>> _coreLossesModels;

            std::vector<double> _magneticFluxDensities = {0.01, 0.025, 0.05, 0.1, 0.2};
            std::vector<double> _frequencies = {20000, 50000, 100000, 250000, 500000};
        public:
            MagneticCoreFilterVolumetricLosses(std::optional<CoreLossesModels> preferredModel = std::nullopt) {
                if (preferredModel) {
                    _coreLossesModels.push_back(std::pair<CoreLossesModels, std::shared_ptr<CoreLossesModel>>{preferredModel.value(), CoreLossesModel::factory(preferredModel.value())});
                }
                // _coreLossesModelNames = {Defaults().coreLossesModelDefault, CoreLossesModels::LOSS_FACTOR, CoreLossesModels::PROPRIETARY, CoreLossesModels::IGSE, CoreLossesModels::ROSHEN};
                _coreLossesModelNames = {Defaults().coreLossesModelDefault, CoreLossesModels::LOSS_FACTOR, CoreLossesModels::PROPRIETARY, CoreLossesModels::IGSE};
                for (auto modelName : _coreLossesModelNames) {
                    _coreLossesModels.push_back(std::pair<CoreLossesModels, std::shared_ptr<CoreLossesModel>>{modelName, CoreLossesModel::factory(modelName)});
                }
            }
            std::vector<std::pair<CoreMaterial, double>> filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, std::map<std::string, std::string> models, double weight=1);
            double calculate_average_volumetric_losses(CoreMaterial coreMaterial, double temperature, std::map<std::string, std::string> models);
    };
    
    class MagneticCoreFilterResistivity : public MagneticCoreFilter {
        public:
            std::vector<std::pair<CoreMaterial, double>> filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, double weight=1);
    };
};


} // namespace OpenMagnetics