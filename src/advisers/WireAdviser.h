#pragma once
#include "constructive_models/Coil.h"
#include "constructive_models/Wire.h"
#include "Defaults.h"
#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {

void from_json(const json & j, WireSolidInsulationRequirements & x);
void to_json(json & j, const WireSolidInsulationRequirements & x);

inline void from_json(const json & j, WireSolidInsulationRequirements& x) {
    x.set_minimum_number_layers(get_stack_optional<int64_t>(j, "minimumNumberLayers"));
    x.set_minimum_grade(get_stack_optional<int64_t>(j, "minimumGrade"));
    x.set_maximum_number_layers(get_stack_optional<int64_t>(j, "maximumNumberLayers"));
    x.set_maximum_grade(get_stack_optional<int64_t>(j, "maximumGrade"));
    x.set_minimum_breakdown_voltage(j.at("minimumBreakdownVoltage").get<double>());
}

inline void to_json(json & j, const WireSolidInsulationRequirements & x) {
    j = json::object();
    j["minimumNumberLayers"] = x.get_minimum_number_layers();
    j["minimumGrade"] = x.get_minimum_grade();
    j["maximumNumberLayers"] = x.get_maximum_number_layers();
    j["maximumGrade"] = x.get_maximum_grade();
    j["minimumBreakdownVoltage"] = x.get_minimum_breakdown_voltage();
}


class WireAdviser {
    protected:
        double _maximumEffectiveCurrentDensity;
        std::optional<WireSolidInsulationRequirements> _wireSolidInsulationRequirements;
        std::optional<WireStandard> _commonWireStandard;
        int _maximumNumberParallels;
        double _maximumOuterAreaProportion;
        std::string _log;
    public:

        WireAdviser(double maximumEffectiveCurrentDensity, double maximumNumberParallels) {
            _maximumEffectiveCurrentDensity = maximumEffectiveCurrentDensity;
            _maximumNumberParallels = maximumNumberParallels;
        }
        WireAdviser() {
            auto defaults = Defaults();
            _maximumEffectiveCurrentDensity = defaults.maximumEffectiveCurrentDensity;
            _maximumNumberParallels = defaults.maximumNumberParallels;
        }
        virtual ~WireAdviser() = default;

        void logEntry(std::string entry);
        
        std::string read_log() {
            return _log;
        }

        void set_maximum_effective_current_density(double maximumEffectiveCurrentDensity) {
            _maximumEffectiveCurrentDensity = maximumEffectiveCurrentDensity;
        }
        void set_wire_solid_insulation_requirements(WireSolidInsulationRequirements wireSolidInsulationRequirements) {
            _wireSolidInsulationRequirements = wireSolidInsulationRequirements;
        }
        void set_maximum_number_parallels(int maximumNumberParallels) {
            _maximumNumberParallels = maximumNumberParallels;
        }
        double get_maximum_area_proportion() {
            return _maximumOuterAreaProportion;
        }
        void set_common_wire_standard(std::optional<WireStandard> commonWireStandard) {
            _commonWireStandard = commonWireStandard;
        }
        std::optional<WireStandard> get_common_wire_standard() {
            return _commonWireStandard;
        }
        std::vector<std::pair<CoilFunctionalDescription, double>> get_advised_wire(CoilFunctionalDescription coilFunctionalDescription,
                                                                                   Section section,
                                                                                   SignalDescriptor current,
                                                                                   double temperature,
                                                                                   uint8_t numberSections,
                                                                                   size_t maximumNumberResults=1);
        std::vector<std::pair<CoilFunctionalDescription, double>> get_advised_wire(std::vector<Wire>* wires,
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

        std::vector<std::pair<CoilFunctionalDescription, double>> filter_by_effective_resistance(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils,
                                                                                                 SignalDescriptor current,
                                                                                                 double temperature);

        std::vector<std::pair<CoilFunctionalDescription, double>> filter_by_proximity_factor(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils,
                                                                                                 SignalDescriptor current,
                                                                                                 double temperature);

        std::vector<std::pair<CoilFunctionalDescription, double>> filter_by_solid_insulation_requirements(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils,
                                                                                                 WireSolidInsulationRequirements wireSolidInsulationRequirements);

        std::vector<std::pair<CoilFunctionalDescription, double>> create_dataset(CoilFunctionalDescription coilFunctionalDescription,
                                                                                 std::vector<Wire>* wires,
                                                                                 Section section,
                                                                                 SignalDescriptor current,
                                                                                 double temperature);
        void expand_wires_dataset_with_parallels(std::vector<CoilFunctionalDescription>* coilFunctionalDescriptions);
        void set_maximum_area_proportion(std::vector<std::pair<CoilFunctionalDescription, double>>* unfilteredCoils, Section section, uint8_t numberSections);
    
};


} // namespace OpenMagnetics