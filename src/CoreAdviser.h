#pragma once
#include "Constants.h"
#include "Defaults.h"
#include "CoilWrapper.h"
#include "OutputsWrapper.h"
#include "InputsWrapper.h"
#include "CoreWrapper.h"
#include "MasWrapper.h"
#include <cmath>
#include <MAS.hpp>


namespace OpenMagnetics {

class CoreAdviser {
    public: 
        enum class CoreAdviserFilters : int {
            AREA_PRODUCT, 
            ENERGY_STORED, 
            COST, 
            EFFICIENCY,
            DIMENSIONS
        };
    protected:
        std::map<std::string, std::string> _models;
        std::string _log;
        bool _uniqueCoreShapes = false;
        std::map<CoreAdviserFilters, double> _weights;
        double _averageMarginInWindingWindow = 0;

        void logEntry(std::string entry) {
            _log += entry + "\n";
        }

    public:

        std::map<CoreAdviserFilters, std::map<std::string, bool>> _filterConfiguration{
                { CoreAdviserFilters::AREA_PRODUCT,          { {"invert", true}, {"log", true} } },
                { CoreAdviserFilters::ENERGY_STORED,         { {"invert", true}, {"log", true} } },
                { CoreAdviserFilters::COST,                  { {"invert", true}, {"log", true} } },
                { CoreAdviserFilters::EFFICIENCY,            { {"invert", true}, {"log", true} } },
                { CoreAdviserFilters::DIMENSIONS,            { {"invert", true}, {"log", true} } },
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
        void set_average_margin_in_winding_window(double value) {
            _averageMarginInWindingWindow = value;
        }
        double get_average_margin_in_winding_window() {
            return _averageMarginInWindingWindow;
        }

        std::vector<std::pair<MasWrapper, double>> get_advised_core(InputsWrapper inputs, size_t maximumNumberResults=1);
        std::vector<std::pair<MasWrapper, double>> get_advised_core(InputsWrapper inputs, std::vector<CoreWrapper>* cores, size_t maximumNumberResults=1);
        std::vector<std::pair<MasWrapper, double>> get_advised_core(InputsWrapper inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumNumberResults=1);
        std::vector<std::pair<MasWrapper, double>> get_advised_core(InputsWrapper inputs, std::map<CoreAdviserFilters, double> weights, std::vector<CoreWrapper>* cores, size_t maximumNumberResults=1);
        std::vector<std::pair<MasWrapper, double>> get_advised_core(InputsWrapper inputs, std::vector<CoreWrapper>* cores, size_t maximumNumberResults, size_t maximumNumberCores);
        std::vector<std::pair<MasWrapper, double>> get_advised_core(InputsWrapper inputs, std::map<CoreAdviserFilters, double> weights, std::vector<CoreWrapper>* cores, size_t maximumNumberResults, size_t maximumNumberCores);

        std::vector<std::pair<MasWrapper, double>> apply_filters(std::vector<std::pair<MasWrapper, double>>* masMagnetics, InputsWrapper inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults);
        std::vector<std::pair<MasWrapper, double>> create_mas_dataset(InputsWrapper inputs, std::vector<CoreWrapper>* cores, bool includeStacks);
        void expand_mas_dataset_with_stacks(InputsWrapper inputs, std::vector<CoreWrapper>* cores, std::vector<std::pair<MasWrapper, double>>* masMagnetics);
    
    class MagneticCoreFilter {
        public:
            std::map<CoreAdviserFilters, std::map<std::string, double>>* _scorings;
            std::map<CoreAdviserFilters, std::map<std::string, bool>>* _validScorings;
            std::map<CoreAdviserFilters, std::map<std::string, bool>>* _filterConfiguration;
            double* _averageMarginInWindingWindow;

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
            void set_average_margin_in_winding_window(double* averageMarginInWindingWindow) {
                _averageMarginInWindingWindow = averageMarginInWindingWindow;
            }
            MagneticCoreFilter(){
            }
            std::vector<std::pair<MasWrapper, double>> filter_magnetics(std::vector<MasWrapper> unfilteredMasMagnetics, InputsWrapper inputs, double weight=1);
    };
    
    class MagneticCoreFilterAreaProduct : public MagneticCoreFilter {
        public:
            std::vector<std::pair<MasWrapper, double>> filter_magnetics(std::vector<std::pair<MasWrapper, double>>* unfilteredMasMagnetics, InputsWrapper inputs, double weight=1, bool firstFilter=false);
    };
    
    class MagneticCoreFilterEnergyStored : public MagneticCoreFilter {
        public:
            std::vector<std::pair<MasWrapper, double>> filter_magnetics(std::vector<std::pair<MasWrapper, double>>* unfilteredMasMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models, double weight=1, bool firstFilter=false);
    };
    
    class MagneticCoreFilterCost : public MagneticCoreFilter {
        public:
            std::vector<std::pair<MasWrapper, double>> filter_magnetics(std::vector<std::pair<MasWrapper, double>>* unfilteredMasMagnetics, InputsWrapper inputs, double weight=1, bool firstFilter=false);
    };
    
    class MagneticCoreFilterLosses : public MagneticCoreFilter {
        public:
            std::vector<std::pair<MasWrapper, double>> filter_magnetics(std::vector<std::pair<MasWrapper, double>>* unfilteredMasMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models, double weight=1, bool firstFilter=false);
    };
    
    class MagneticCoreFilterDimensions : public MagneticCoreFilter {
        public:
            std::vector<std::pair<MasWrapper, double>> filter_magnetics(std::vector<std::pair<MasWrapper, double>>* unfilteredMasMagnetics, double weight=1, bool firstFilter=false);
    };

};


} // namespace OpenMagnetics