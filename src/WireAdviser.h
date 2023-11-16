#pragma once
#include "WireWrapper.h"
#include "Defaults.h"
#include <MAS.hpp>


namespace OpenMagnetics {

class WireAdviser {
    protected:
        bool _includeFoil;
        bool _includeRectangular;
        bool _includeLitz;
        bool _includeRound;
        double _maximumEffectiveCurrentDensity;
        int _maximumNumberParallels;
    public:
        WireAdviser(double maximumEffectiveCurrentDensity, double maximumNumberParallels, bool includeRectangular=true, bool includeFoil=false) {
            _includeFoil = includeFoil;
            _includeRectangular = includeRectangular;
            _maximumEffectiveCurrentDensity = maximumEffectiveCurrentDensity;
            _maximumNumberParallels = maximumNumberParallels;
        }
        WireAdviser(bool includeRectangular=true, bool includeFoil=false, bool includeLitz=true, bool includeRound=true) {
            auto defaults = Defaults();
            _includeFoil = includeFoil;
            _includeRectangular = includeRectangular;
            _includeLitz = includeLitz;
            _includeRound = includeRound;
            _maximumEffectiveCurrentDensity = defaults.maximumEffectiveCurrentDensity;
            _maximumNumberParallels = defaults.maximumNumberParallels;
        }
        virtual ~WireAdviser() = default;

        void set_maximum_effective_current_density(double maximumEffectiveCurrentDensity) {
            _maximumEffectiveCurrentDensity = maximumEffectiveCurrentDensity;
        }
        void set_maximum_number_parallels(int maximumNumberParallels) {
            _maximumNumberParallels = maximumNumberParallels;
        }
        void enableFoil(bool enable) {
            _includeFoil = enable;
        }
        void enableRectangular(bool enable) {
            _includeRectangular = enable;
        }
        void enableLitz(bool enable) {
            _includeLitz = enable;
        }
        void enableRound(bool enable) {
            _includeRound = enable;
        }
        std::vector<std::pair<CoilFunctionalDescription, double>> get_advised_wire(CoilFunctionalDescription coilFunctionalDescription,
                                                                                   Section section,
                                                                                   SignalDescriptor current,
                                                                                   double temperature,
                                                                                   uint8_t numberSections,
                                                                                   size_t maximumNumberResults=1);
        std::vector<std::pair<CoilFunctionalDescription, double>> get_advised_wire(std::vector<WireWrapper>* wires,
                                                                                   CoilFunctionalDescription coilFunctionalDescription,
                                                                                   Section section,
                                                                                   SignalDescriptor current,
                                                                                   double temperature,
                                                                                   uint8_t numberSections,
                                                                                   size_t maximumNumberResults=1);

        std::vector<std::pair<CoilFunctionalDescription, double>> filter_by_area_no_parallels(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils,
                                                                                              Section section);

        std::vector<std::pair<CoilFunctionalDescription, double>> filter_by_area_with_parallels(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils,
                                                                                                Section section,
                                                                                                double numberSections);

        std::vector<std::pair<CoilFunctionalDescription, double>> filter_by_skin_depth_resistance(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils,
                                                                                                  SignalDescriptor current,
                                                                                                  double temperature);

        std::vector<std::pair<CoilFunctionalDescription, double>> filter_by_effective_resistance(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils,
                                                                                                 SignalDescriptor current,
                                                                                                 double temperature);

        std::vector<std::pair<CoilFunctionalDescription, double>> filter_by_proximity_factor(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils,
                                                                                                 SignalDescriptor current,
                                                                                                 double temperature);

        std::vector<std::pair<CoilFunctionalDescription, double>> create_dataset(CoilFunctionalDescription coilFunctionalDescription,
                                                                                 std::vector<WireWrapper>* wires,
                                                                                 Section section,
                                                                                 SignalDescriptor current,
                                                                                 double temperature,
                                                                                 bool includeExtraParallels);
        void expand_wires_dataset_with_parallels(std::vector<CoilFunctionalDescription>* coilFunctionalDescriptions);
    
};


} // namespace OpenMagnetics