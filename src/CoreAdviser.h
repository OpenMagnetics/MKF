#pragma once
#include "Constants.h"
#include "Defaults.h"
#include "CoilWrapper.h"
#include "OutputsWrapper.h"
#include "InputsWrapper.h"
#include "CoreWrapper.h"
#include "MasWrapper.h"
#include <MAS.hpp>


namespace OpenMagnetics {

class CoreAdviser {
    public: 
        enum class CoreAdviserFilters : int {AREA_PRODUCT, 
            ENERGY_STORED, 
            WINDING_WINDOW_AREA, 
            CORE_LOSSES, 
            CORE_TEMPERATURE, 
            DIMENSIONS
        };
    private:
        std::map<std::string, std::string> _models;
        bool _includeToroids;
        std::string _log;
        std::map<CoreAdviserFilters, double> _weights;

        void log(std::string entry) {
            _log += entry + "\n";
        }
    public:
        std::map<CoreAdviserFilters, std::map<std::string, double>> _scorings;
        CoreAdviser(std::map<std::string, std::string> models, bool includeToroids=true) {
            auto defaults = OpenMagnetics::Defaults();
            _includeToroids = includeToroids;
            _models = models;
            if (models.find("gapReluctance") == models.end()) {
                _models["gapReluctance"] = magic_enum::enum_name(defaults.reluctanceModelDefault);
            }
        }
        CoreAdviser(bool includeToroids=true) {
            auto defaults = OpenMagnetics::Defaults();
            _includeToroids = includeToroids;
            _models["gapReluctance"] = magic_enum::enum_name(defaults.reluctanceModelDefault);
            _models["coreLosses"] = magic_enum::enum_name(defaults.coreLossesModelDefault);
            _models["coreTemperature"] = magic_enum::enum_name(defaults.coreTemperatureModelDefault);
        }
        std::string read_log() {
            return _log;
        }
        std::map<std::string, std::map<CoreAdviserFilters, double>> get_scorings(bool weighted = false) {
            std::map<std::string, std::map<CoreAdviserFilters, double>> invertedScorings;
            for (auto const& [filter, aux] : _scorings) {
                double maximumScoring = (*std::max_element(aux.begin(), aux.end(),
                                             [](const std::pair<std::string, double> &p1,
                                                const std::pair<std::string, double> &p2)
                                             {
                                                 return p1.second < p2.second;
                                             })).second; 
                double minimumScoring = (*std::min_element(aux.begin(), aux.end(),
                                             [](const std::pair<std::string, double> &p1,
                                                const std::pair<std::string, double> &p2)
                                             {
                                                 return p1.second < p2.second;
                                             })).second; 

                for (auto const& [name, scoring] : aux) {
                    if (weighted)
                        invertedScorings[name][filter] = _weights[filter] * (scoring - minimumScoring) / (maximumScoring - minimumScoring);
                    else
                        invertedScorings[name][filter] = (scoring - minimumScoring) / (maximumScoring - minimumScoring);
                }
            }
            return invertedScorings;
        }

        std::vector<MasWrapper> get_advised_core(InputsWrapper inputs, size_t maximumNumberResults=1);
        std::vector<MasWrapper> get_advised_core(InputsWrapper inputs, std::vector<CoreWrapper> cores, size_t maximumNumberResults=1);
        std::vector<MasWrapper> get_advised_core(InputsWrapper inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumNumberResults=1);
        std::vector<MasWrapper> get_advised_core(InputsWrapper inputs, std::map<CoreAdviserFilters, double> weights, std::vector<CoreWrapper> cores, size_t maximumNumberResults=1);
        std::vector<MasWrapper> apply_filters(std::vector<std::pair<MasWrapper, double>> masMagnetics, InputsWrapper inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults);
    
    class MagneticCoreFilter {
        std::map<CoreAdviserFilters, std::map<std::string, double>>* _scorings;

        public:

            void add_scoring(std::string name, CoreAdviser::CoreAdviserFilters filter, double scoring) {
                // if (!(*_scorings).contains(name)) {
                    // (*_scorings)[name] = std::map<CoreAdviserFilters, double>{filter,;
                // }
                // std::cout << scoring << std::endl;
                (*_scorings)[filter][name] = scoring;
                // std::cout << (*_scorings)[name][filter] << std::endl;

            }
            MagneticCoreFilter(std::map<CoreAdviserFilters, std::map<std::string, double>>* scorings){
                _scorings = scorings;
            }
            std::vector<std::pair<MasWrapper, double>> filter_magnetics(std::vector<MasWrapper> unfilteredMasMagnetics, InputsWrapper inputs, double weight=1);
    };
    
    class MagneticCoreFilterAreaProduct : public MagneticCoreFilter {
        public:
            std::vector<std::pair<MasWrapper, double>> filter_magnetics(std::vector<std::pair<MasWrapper, double>> unfilteredMasMagnetics, InputsWrapper inputs, double weight=1);
    };
    
    class MagneticCoreFilterEnergyStored : public MagneticCoreFilter {
        public:
            std::vector<std::pair<MasWrapper, double>> filter_magnetics(std::vector<std::pair<MasWrapper, double>> unfilteredMasMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models, double weight=1);
    };
    
    class MagneticCoreFilterWindingWindowArea : public MagneticCoreFilter {
        public:
            std::vector<std::pair<MasWrapper, double>> filter_magnetics(std::vector<std::pair<MasWrapper, double>> unfilteredMasMagnetics, InputsWrapper inputs, double weight=1);
    };
    
    class MagneticCoreFilterCoreLosses : public MagneticCoreFilter {
        public:
            std::vector<std::pair<MasWrapper, double>> filter_magnetics(std::vector<std::pair<MasWrapper, double>> unfilteredMasMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models, double weight=1);
    };
    
    class MagneticCoreFilterCoreTemperature : public MagneticCoreFilter {
        public:
            std::vector<std::pair<MasWrapper, double>> filter_magnetics(std::vector<std::pair<MasWrapper, double>> unfilteredMasMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models, double weight=1);
    };
    
    class MagneticCoreFilterDimensions : public MagneticCoreFilter {
        public:
            std::vector<std::pair<MasWrapper, double>> filter_magnetics(std::vector<std::pair<MasWrapper, double>> unfilteredMasMagnetics, double weight=1);
    };

};


} // namespace OpenMagnetics