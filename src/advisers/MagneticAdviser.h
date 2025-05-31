#pragma once
#include "support/Utils.h"
#include "advisers/MagneticFilter.h"
#include "advisers/CoreAdviser.h"
#include "advisers/CoilAdviser.h"
#include "constructive_models/Mas.h"
#include "Definitions.h"
#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {

class MagneticAdviser{
    public:

        std::map<MagneticFilters, std::map<std::string, double>> _scorings;
        std::map<MagneticFilters, std::shared_ptr<MagneticFilter>> _filters;
        std::vector<MagneticFilterOperation> _loadedFilterFlow;
        std::vector<MagneticFilterOperation> _defaultCustomMagneticFilterFlow{
            MagneticFilterOperation(MagneticFilters::COST, true, true, 1.0),
            MagneticFilterOperation(MagneticFilters::LOSSES, true, true, 1.0),
            MagneticFilterOperation(MagneticFilters::DIMENSIONS, true, true, 1.0),
        };

        std::vector<MagneticFilterOperation> _defaultCatalogMagneticFilterFlow{
            MagneticFilterOperation(MagneticFilters::TURNS_RATIOS, true, false, true, 1.0),
            MagneticFilterOperation(MagneticFilters::MAXIMUM_DIMENSIONS, true, false, 1.0),
            MagneticFilterOperation(MagneticFilters::SATURATION, true, false, 1.0),
            MagneticFilterOperation(MagneticFilters::DC_CURRENT_DENSITY, true, false, 1.0),
            MagneticFilterOperation(MagneticFilters::EFFECTIVE_CURRENT_DENSITY, true, false, 1.0),
            MagneticFilterOperation(MagneticFilters::IMPEDANCE, true, false, 1.0),
            MagneticFilterOperation(MagneticFilters::MAGNETIZING_INDUCTANCE, true, false, 1.0),
        };
        bool _simulateResults = true;

        MagneticAdviser() {
        }
        MagneticAdviser(bool simulateResults) {
            _simulateResults = simulateResults;
        }

        void add_scoring(std::string name, MagneticFilters filter, double scoring) {
            if (scoring != -1) {
                _scorings[filter][name] = scoring;
            }
        }

        std::vector<std::pair<Mas, double>> get_advised_magnetic(Inputs inputs, size_t maximumNumberResults=1);
        std::vector<std::pair<Mas, double>> get_advised_magnetic(Inputs inputs, std::map<MagneticFilters, double> weights, size_t maximumNumberResults);
        std::vector<std::pair<Mas, double>> get_advised_magnetic(Inputs inputs, std::vector<MagneticFilterOperation> filterFlow, size_t maximumNumberResults);
        std::vector<std::pair<Mas, double>> get_advised_magnetic(Inputs inputs, std::vector<Magnetic> catalogMagnetics, size_t maximumNumberResults=1, bool strict=true);
        std::vector<std::pair<Mas, double>> get_advised_magnetic(Inputs inputs, std::vector<Magnetic> catalogMagnetics, std::vector<MagneticFilterOperation> filterFlow, size_t maximumNumberResults=1, bool strict=true);
        std::vector<std::pair<Mas, double>> get_advised_magnetic(std::vector<Mas> catalogMagneticsWithInputs, std::vector<MagneticFilterOperation> filterFlow, size_t maximumNumberResults=1, bool strict=true);
        std::vector<std::pair<Mas, double>> score_magnetics(std::vector<Mas> masMagneticsWithCoil, std::vector<MagneticFilterOperation> filterFlow);
        void normalize_scoring(std::vector<std::pair<Mas, double>>* masMagneticsWithScoring, std::vector<double> scoring, double weight, std::map<std::string, bool> filterConfiguration);
        void normalize_scoring(std::vector<std::pair<Mas, double>>* masMagneticsWithScoring, std::vector<double> scoring, MagneticFilterOperation filterConfiguration);
        static void preview_magnetic(Mas mas);
        void load_filter_flow(std::vector<MagneticFilterOperation> flow);
        std::map<std::string, std::map<MagneticFilters, double>> get_scorings();
};


} // namespace OpenMagnetics