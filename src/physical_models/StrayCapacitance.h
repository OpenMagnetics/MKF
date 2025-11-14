#pragma once
#include "Defaults.h"
#include "constructive_models/Magnetic.h"
#include "support/Utils.h"
#include <MAS.hpp>
#include "Models.h"

using namespace MAS;

namespace OpenMagnetics {


class StrayCapacitanceModel {
    private:
    protected:

    public:
        std::string methodName = "Default";
        static std::shared_ptr<StrayCapacitanceModel> factory(StrayCapacitanceModels modelName);
        std::vector<double> preprocess_data(Turn firstTurn, Wire firstWire, Turn secondTurn, Wire secondWire, Coil coil);
        virtual double calculate_static_capacitance_between_two_turns(double insulationThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double epsilonD, double epsilonF) = 0;
};

// Based on "Self-Capacitance of Inductors" by Antonio Massarini
// https://sci-hub.st/https://ieeexplore.ieee.org/document/602562
class StrayCapacitanceMassariniModel : public StrayCapacitanceModel {
    public:
        std::string methodName = "Massarini";
        double calculate_static_capacitance_between_two_turns(double insulationThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double epsilonD, double epsilonF);

};

// Based on "Equivalent capacitances of transformer windings" by W. T. Duerdoth
class StrayCapacitanceDuerdothModel : public StrayCapacitanceModel {
    public:
        std::string methodName = "Duerdoth";
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
            _model = StrayCapacitanceModel::factory(StrayCapacitanceModels::ALBACH);
        };
        virtual ~StrayCapacitance() = default;


        static std::vector<std::pair<Turn, size_t>> get_surrounding_turns(Turn currentTurn, std::vector<Turn> turnsDescription);
        static StrayCapacitanceOutput calculate_voltages_per_turn(Coil coil, OperatingPoint operatingPoint);
        static StrayCapacitanceOutput calculate_voltages_per_turn(Coil coil, std::map<std::string, double> voltageRmsPerWinding);
        static std::vector<Layer> get_insulation_layers_between_two_turns(Turn firstTurn, Turn secondTurn, Coil coil);
        double calculate_static_capacitance_between_two_turns(Turn firstTurn, Wire firstWire, Turn secondTurn, Wire secondWire, Coil coil);

        std::map<std::pair<size_t, size_t>, double> calculate_capacitance_among_turns(Coil coil);

        std::map<std::pair<std::string, std::string>, double> calculate_capacitance_among_windings(Coil coil);
        static std::map<std::string, double> calculate_capacitance_matrix(double energy, double voltageDrop, double relativeTurnsRatio);
        std::map<std::pair<std::string, std::string>, double> calculate_maxwell_capacitance_matrix(Coil coil);

};



class StrayCapacitanceOneLayer{

    public:

        StrayCapacitanceOneLayer(){
        };
        virtual ~StrayCapacitanceOneLayer() = default;
        double calculate_capacitance(Coil coil);


};
} // namespace OpenMagnetics
