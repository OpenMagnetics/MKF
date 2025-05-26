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
        enum class MagneticAdviserFilters : int {
            COST, 
            EFFICIENCY,
            DIMENSIONS
        };
    private:
        std::map<MagneticAdviserFilters, double> _weights;
    public:

        std::map<MagneticAdviserFilters, std::map<std::string, bool>> _filterConfiguration{
                { MagneticAdviserFilters::COST,                  { {"invert", true}, {"log", false} } },
                { MagneticAdviserFilters::EFFICIENCY,            { {"invert", true}, {"log", false} } },
                { MagneticAdviserFilters::DIMENSIONS,            { {"invert", true}, {"log", false} } },
            };
        std::map<MagneticAdviserFilters, std::map<std::string, double>> _scorings;
        std::map<MagneticFilters, std::shared_ptr<MagneticFilter>> _filters;
        std::vector<MagneticFilterOperation> _defaultCustomMagneticFilterFlow{
            MagneticFilterOperation(MagneticFilters::COST, true, true, 1.0),
            MagneticFilterOperation(MagneticFilters::LOSSES, true, true, 1.0),
            MagneticFilterOperation(MagneticFilters::DIMENSIONS, true, true, 1.0),
        };

        MagneticAdviser() {
            for (auto filterConfiguration : _defaultCustomMagneticFilterFlow) {
                MagneticFilters filterEnum = filterConfiguration.get_filter();
                std::shared_ptr<MagneticFilter> filter = MagneticFilter::factory(filterEnum);
                _filters[filterEnum] = filter;
            }
        }

        void add_scoring(std::string name, MagneticAdviser::MagneticAdviserFilters filter, double scoring) {
            if (scoring != -1) {
                _scorings[filter][name] = scoring;
            }
        }

        std::vector<std::pair<Mas, double>> get_advised_magnetic(Inputs inputs, size_t maximumNumberResults=1);
        std::vector<std::pair<Mas, double>> get_advised_magnetic(Inputs inputs, std::map<MagneticAdviserFilters, double> weights, size_t maximumNumberResults);
        std::vector<std::pair<Mas, double>> get_advised_magnetic(Inputs inputs, std::vector<Magnetic> catalogMagnetics, size_t maximumNumberResults=1);
        std::vector<std::pair<Mas, double>> score_magnetics(std::vector<Mas> masMagneticsWithCoil, std::map<MagneticAdviserFilters, double> weights);
        void normalize_scoring(std::vector<std::pair<Mas, double>>* masMagneticsWithScoring, std::vector<double> scoring, double weight, std::map<std::string, bool> filterConfiguration);
        static void preview_magnetic(Mas mas);

        std::map<std::string, std::map<MagneticAdviserFilters, double>> get_scorings(){
            return get_scorings(false);
        }

        std::map<std::string, std::map<MagneticAdviserFilters, double>> get_scorings(bool weighted);
};


} // namespace OpenMagnetics