#pragma once
#include "Constants.h"
#include "Defaults.h"
#include "WindingWrapper.h"
#include "InputsWrapper.h"
#include "CoreWrapper.h"
#include <MAS.hpp>


namespace OpenMagnetics {

class CoreAdviser {
  private:
    std::map<std::string, std::string> _models;
  protected:
  public:
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
    CoreWrapper get_advised_core(InputsWrapper inputs);

    class MagneticCoreFilter {
      public:
        static std::pair<std::vector<Magnetic>, std::map<std::string, double>> filter_magnetics(std::vector<Magnetic> unfilteredMagnetics, InputsWrapper inputs);
    };

    class MagneticCoreFilterAreaProduct : public MagneticCoreFilter {
      public:
        static std::pair<std::vector<Magnetic>, std::map<std::string, double>> filter_magnetics(std::vector<Magnetic> unfilteredMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models);
    };

    class MagneticCoreFilterEnergyStored : public MagneticCoreFilter {
      public:
        static std::pair<std::vector<Magnetic>, std::map<std::string, double>> filter_magnetics(std::vector<Magnetic> unfilteredMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models);
    };

    class MagneticCoreFilterWindingWindowArea : public MagneticCoreFilter {
      public:
        static std::pair<std::vector<Magnetic>, std::map<std::string, double>> filter_magnetics(std::vector<Magnetic> unfilteredMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models);
    };

    class MagneticCoreFilterCoreLosses : public MagneticCoreFilter {
      public:
        static std::pair<std::vector<Magnetic>, std::map<std::string, double>> filter_magnetics(std::vector<Magnetic> unfilteredMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models);
    };

    class MagneticCoreFilterCoreTemperature : public MagneticCoreFilter {
      public:
        static std::pair<std::vector<Magnetic>, std::map<std::string, double>> filter_magnetics(std::vector<Magnetic> unfilteredMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models, std::map<std::string, double> precalculatedCoreLosses);
    };

};


} // namespace OpenMagnetics