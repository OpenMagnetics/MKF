#pragma once
#include "Defaults.h"
#include "MagneticWrapper.h"
#include "Utils.h"
#include <MAS.hpp>
#include "Models.h"

namespace OpenMagnetics {


class StrayCapacitanceModel {
    private:
    protected:

    public:
        std::string methodName = "Default";
        static std::shared_ptr<StrayCapacitanceModel> factory(StrayCapacitanceModels modelName);
        std::vector<double> preprocess_data(Turn firstTurn, WireWrapper firstWire, Turn secondTurn, WireWrapper secondWire, CoilWrapper coil);
        virtual double calculate_static_capacitance_between_two_turns(double insulationThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double epsilonD, double epsilonF) = 0;
};

// Based on "Induktivitäten in der Leistungselektronik", pages 49-50, by Manfred Albach
class StrayCapacitanceDuerdothModel : public StrayCapacitanceModel {
    public:
        std::string methodName = "Albach";
        double calculate_static_capacitance_between_two_turns(double insulationThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double epsilonD, double epsilonF);

};

// Based on "Induktivitäten in der Leistungselektronik", pages 49-50, by Manfred Albach
class StrayCapacitanceAlbachModel : public StrayCapacitanceModel {
    public:
        std::string methodName = "Albach";
        double calculate_static_capacitance_between_two_turns(double insulationThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double epsilonD, double epsilonF);

};

// Based on "“Berechnung der kapazitat von spulen, insbesondere in schalenkernen" by K. Koch
// Reproduced in "Using Transformer Parasitics for Resonant Converters—A Review of the Calculation of the Stray Capacitance of Transformers" by Juergen Biela and Johann W. Kolar  
// https://www.pes-publications.ee.ethz.ch/uploads/tx_ethpublications/biela_IEEETrans_ReviewStrayCap.pdf
class StrayCapacitanceKochModel : public StrayCapacitanceModel {
    public:
        std::string methodName = "Koch";
        double calculate_static_capacitance_between_two_turns(double insulationThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double epsilonD, double epsilonF);

};


class StrayCapacitance{
    std::shared_ptr<StrayCapacitanceModel> _model;
    public:

        StrayCapacitance(){
            _model = StrayCapacitanceModel::factory(StrayCapacitanceModels::KOCH);
        };
        virtual ~StrayCapacitance() = default;


        static std::vector<Turn> get_surrounding_turns(Turn currentTurn, std::vector<Turn> turnsDescription);
        static StrayCapacitanceOutput calculate_voltages_per_turn(CoilWrapper coil, OperatingPoint operatingPoint);
        static std::vector<Layer> get_insulation_layers_between_two_turns(Turn firstTurn, Turn secondTurn, CoilWrapper coil);
        double calculate_static_capacitance_between_two_turns(Turn firstTurn, WireWrapper firstWire, Turn secondTurn, WireWrapper secondWire, CoilWrapper coil);

        // This function, if scaled properly, should give the 6 capacitors equivalent, following Kolar's paper for connecting layers, but with turns
        std::map<std::pair<std::string, std::string>, double> calculate_capacitance_among_turns(CoilWrapper coil);
        std::map<std::pair<std::string, std::string>, double> calculate_capacitance_among_layers(CoilWrapper coil);

        // This function, if done properly, should give the 3 capacitors equivalent
        std::map<std::pair<std::string, std::string>, double> calculate_capacitance_among_windings(CoilWrapper coil);

};
} // namespace OpenMagnetics
