#pragma once
#include "WireWrapper.h"
#include "Defaults.h"
#include <MAS.hpp>


namespace OpenMagnetics {

class WireSolidInsulationRequirements {
    public:
    WireSolidInsulationRequirements() = default;
    virtual ~WireSolidInsulationRequirements() = default;

    private:
    std::optional<int64_t> minimumNumberLayers;
    std::optional<int64_t> minimumGrade;
    double minimumBreakdownVoltage = 0;

    public:

    std::optional<int64_t> get_minimum_number_layers() const { return minimumNumberLayers; }
    void set_minimum_number_layers(const int64_t & value) { this->minimumNumberLayers = value; }

    std::optional<int64_t> get_minimum_grade() const { return minimumGrade; }
    void set_minimum_grade(const int64_t & value) { this->minimumGrade = value; }

    double get_minimum_breakdown_voltage() const { return minimumBreakdownVoltage; }
    void set_minimum_breakdown_voltage(const double & value) { this->minimumBreakdownVoltage = value; }
};

class WireAdviser {
    protected:
        bool _includeFoil;
        bool _includeRectangular;
        bool _includeLitz;
        bool _includeRound;
        double _maximumEffectiveCurrentDensity;
        std::optional<WireSolidInsulationRequirements> _wireSolidInsulationRequirements;
        int _maximumNumberParallels;
        double _maximumOuterAreaProportion;
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
        void set_wire_solid_insulation_requirements(WireSolidInsulationRequirements wireSolidInsulationRequirements) {
            _wireSolidInsulationRequirements = wireSolidInsulationRequirements;
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
        double get_maximum_area_proportion() {
            return _maximumOuterAreaProportion;
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
                                                                                                double numberSections,
                                                                                                bool allowNotFit);

        std::vector<std::pair<CoilFunctionalDescription, double>> filter_by_skin_depth_resistance(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils,
                                                                                                  SignalDescriptor current,
                                                                                                  double temperature);

        std::vector<std::pair<CoilFunctionalDescription, double>> filter_by_effective_resistance(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils,
                                                                                                 SignalDescriptor current,
                                                                                                 double temperature);

        std::vector<std::pair<CoilFunctionalDescription, double>> filter_by_proximity_factor(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils,
                                                                                                 SignalDescriptor current,
                                                                                                 double temperature);

        std::vector<std::pair<CoilFunctionalDescription, double>> filter_by_solid_insulation_requirements(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils,
                                                                                                 WireSolidInsulationRequirements wireSolidInsulationRequirements);

        std::vector<std::pair<CoilFunctionalDescription, double>> create_dataset(CoilFunctionalDescription coilFunctionalDescription,
                                                                                 std::vector<WireWrapper>* wires,
                                                                                 Section section,
                                                                                 SignalDescriptor current,
                                                                                 double temperature);
        void expand_wires_dataset_with_parallels(std::vector<CoilFunctionalDescription>* coilFunctionalDescriptions);
        void set_maximum_area_proportion(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils, Section section, uint8_t numberSections);
    
};


} // namespace OpenMagnetics