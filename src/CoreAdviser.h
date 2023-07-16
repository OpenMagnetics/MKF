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
    private:
        std::map<std::string, std::string> _models;
        bool _includeToroids;
        std::string _log;
    protected:
    public: 
        enum class CoreAdviserFilters : int {AREA_PRODUCT, 
            ENERGY_STORED, 
            WINDING_WINDOW_AREA, 
            CORE_LOSSES, 
            CORE_TEMPERATURE, 
            DIMENSIONS
        };
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

        void log(std::string entry) {
            _log += entry + "\n";
        }
        std::string read_log() {
            return _log;
        }

        std::vector<MasWrapper> get_advised_core(InputsWrapper inputs, size_t maximumNumberResults=1);
        std::vector<MasWrapper> get_advised_core(InputsWrapper inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumNumberResults=1);
        std::vector<MasWrapper> apply_filters(std::vector<std::pair<MasWrapper, double>> masMagnetics, InputsWrapper inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumMagneticsAfterFiltering, size_t maximumNumberResults);
    
    class MagneticCoreFilter {
        public:
            MagneticCoreFilter() {
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