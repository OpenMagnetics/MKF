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
        std::vector<double> preprocess_data_for_round_wires(Turn firstTurn, Wire firstWire, Turn secondTurn, Wire secondWire, std::optional<Coil> coil = std::nullopt);
        virtual double calculate_static_capacitance_between_two_turns(double wireCoatingThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double relativePermittivityWireCoating, double relativePermittivityInsulationLayers) = 0;
};


class StrayCapacitanceParallelPlateModel {
    public:
        std::string methodName = "ParallelPlate";
        double calculate_static_capacitance_between_two_turns(double overlappingDimension, double averageTurnLength, double distanceThroughLayers, double relativePermittivityInsulationLayers);
        std::vector<double> preprocess_data_for_planar_wires(Turn firstTurn, Wire firstWire, Turn secondTurn, Wire secondWire);

};

// Based on "Self-Capacitance of Inductors" by Antonio Massarini
// https://sci-hub.st/https://ieeexplore.ieee.org/document/602562
class StrayCapacitanceMassariniModel : public StrayCapacitanceModel {
    public:
        // Sets the INHERITED methodName rather than redeclaring one. A redeclaration here shadows
        // the base member, and StrayCapacitance holds the model through a
        // shared_ptr<StrayCapacitanceModel>, so `_model->methodName` bound to the BASE member and
        // every result reported "Default" -- including results where the requested model really
        // had been applied. See ABT #950.
        StrayCapacitanceMassariniModel() { methodName = "Massarini"; }
        double calculate_static_capacitance_between_two_turns(double wireCoatingThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double relativePermittivityWireCoating, double relativePermittivityInsulationLayers);

};

// Based on "Equivalent capacitances of transformer windings" by W. T. Duerdoth
class StrayCapacitanceDuerdothModel : public StrayCapacitanceModel {
    public:
        // Sets the INHERITED methodName rather than redeclaring one. A redeclaration here shadows
        // the base member, and StrayCapacitance holds the model through a
        // shared_ptr<StrayCapacitanceModel>, so `_model->methodName` bound to the BASE member and
        // every result reported "Default" -- including results where the requested model really
        // had been applied. See ABT #950.
        StrayCapacitanceDuerdothModel() { methodName = "Duerdoth"; }
        double calculate_static_capacitance_between_two_turns(double wireCoatingThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double relativePermittivityWireCoating, double relativePermittivityInsulationLayers);

};

// Based on "Induktivitäten in der Leistungselektronik", pages 49-50, by Manfred Albach
class StrayCapacitanceAlbachModel : public StrayCapacitanceModel {
    public:
        // Sets the INHERITED methodName rather than redeclaring one. A redeclaration here shadows
        // the base member, and StrayCapacitance holds the model through a
        // shared_ptr<StrayCapacitanceModel>, so `_model->methodName` bound to the BASE member and
        // every result reported "Default" -- including results where the requested model really
        // had been applied. See ABT #950.
        StrayCapacitanceAlbachModel() { methodName = "Albach"; }
        double calculate_static_capacitance_between_two_turns(double wireCoatingThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double relativePermittivityWireCoating, double relativePermittivityInsulationLayers);

};

// Based on "“Berechnung der kapazitat von spulen, insbesondere in schalenkernen" by K. Koch
// Reproduced in "Using Transformer Parasitics for Resonant Converters—A Review of the Calculation of the Stray Capacitance of Transformers" by Juergen Biela and Johann W. Kolar  
// https://www.pes-publications.ee.ethz.ch/uploads/tx_ethpublications/biela_IEEETrans_ReviewStrayCap.pdf
class StrayCapacitanceKochModel : public StrayCapacitanceModel {
    public:
        // Sets the INHERITED methodName rather than redeclaring one. A redeclaration here shadows
        // the base member, and StrayCapacitance holds the model through a
        // shared_ptr<StrayCapacitanceModel>, so `_model->methodName` bound to the BASE member and
        // every result reported "Default" -- including results where the requested model really
        // had been applied. See ABT #950.
        StrayCapacitanceKochModel() { methodName = "Koch"; }
        double calculate_static_capacitance_between_two_turns(double wireCoatingThickness, double averageTurnLength, double conductingRadius, double distanceThroughLayers, double distanceThroughAir, double relativePermittivityWireCoating, double relativePermittivityInsulationLayers);

};


class StrayCapacitance{
    private:
        std::shared_ptr<StrayCapacitanceModel> _model;
        StrayCapacitanceModels _modelName;
        static double calculate_area_between_two_turns_using_diagonals(Turn firstTurn, Turn secondTurn);
        static double calculate_area_between_two_turns_using_vecticals_and_horizontals(Turn firstTurn, Turn secondTurn);
        StrayCapacitanceOutput calculate_capacitance_with_voltages(Coil coil, std::map<std::string, double> voltageRmsPerWinding, std::optional<Core> core = std::nullopt, std::optional<double> frequency = std::nullopt);
    public:

        StrayCapacitance(StrayCapacitanceModels strayCapacitanceModel = StrayCapacitanceModels::ALBACH){
            _model = StrayCapacitanceModel::factory(strayCapacitanceModel);
            _modelName = strayCapacitanceModel;
        };
        virtual ~StrayCapacitance() = default;


        static std::vector<std::pair<Turn, size_t>> get_surrounding_turns(Turn currentTurn, std::vector<Turn> turnsDescription, double globalMinimumGap = -1.0);
        static StrayCapacitanceOutput calculate_voltages_per_turn(Coil coil, OperatingPoint operatingPoint);
        static StrayCapacitanceOutput calculate_voltages_per_turn(Coil coil, std::map<std::string, double> voltageRmsPerWinding);
        static std::vector<Layer> get_insulation_layers_between_two_turns(Turn firstTurn, Turn secondTurn, Coil coil);
        double calculate_static_capacitance_between_two_turns(Turn firstTurn, Wire firstWire, Turn secondTurn, Wire secondWire, std::optional<Coil> coil = std::nullopt);
        double calculate_energy_between_two_turns(Turn firstTurn, Wire firstWire, Turn secondTurn, Wire secondWire, double voltageDrop, std::optional<Coil> coil = std::nullopt);
        double calculate_energy_density_between_two_turns(Turn firstTurn, Wire firstWire, Turn secondTurn, Wire secondWire, double voltageDrop, std::optional<Coil> coil = std::nullopt);
        static double calculate_area_between_two_turns(Turn firstTurn, Turn secondTurn);

        // Capacitance of a single turn to the (equipotential) ferrite core surface,
        // through the dielectric stack wire-enamel | air | core-coating. Reuses the
        // Massarini turn-to-turn formula (inherently bounded — it depends on ln(D0/Dc)
        // set by the wire coating, not an acosh that diverges at contact), with the core
        // coating as the inter-electrode insulation layer, x2 for the wire-to-plane
        // (image) geometry. The core coating gives the finite floor at zero air gap.
        // Building block for the through-core inter-winding path (separated-winding CMC
        // differential mode): summed over a winding's turns by calculate_winding_to_core_capacitance.
        static double calculate_turn_to_core_capacitance(double conductingRadius, double turnLength,
                                                         double wireCoatingThickness, double wireCoatingRelativePermittivity,
                                                         double airGapToCore,
                                                         double coreCoatingThickness, double coreCoatingRelativePermittivity);

        // Total capacitance from one winding to the (equipotential) ferrite core: the
        // parallel sum of its turns' turn-to-core elements. Two of these in series through
        // the core node give the inter-winding capacitance for separated windings.
        static double calculate_winding_to_core_capacitance(Coil coil, Core core, std::string windingName, std::optional<double> frequency = std::nullopt);

        // Inter-winding capacitance between two SEPARATED windings through the floating,
        // equipotential ferrite core (turn -> core -> turn). Energy method: weights each
        // turn-to-core element by the actual per-turn potential (voltagesPerTurn) and
        // solves the floating-core node, per the CPSS 2025 core-potential method — so the
        // turns' potential distribution is accounted for, not a naive parallel sum (which
        // overestimates by ~3x). The two windings carry opposing DM currents, so the
        // second winding's per-turn potentials enter with the opposite sign.
        static double calculate_through_core_capacitance(Coil coil, Core core,
                                                         const std::string& firstWindingName,
                                                         const std::string& secondWindingName,
                                                         const std::vector<double>& voltagesPerTurn,
                                                         std::optional<double> frequency = std::nullopt);

        // Energy stored in ONE winding's turn-to-core elements against the floating core
        // (ABT #848): same per-turn elements and charge-balanced core node as
        // calculate_through_core_capacitance, but for a single winding driven alone —
        // the missing HALF of a winding's self-capacitance. The turn-to-turn chain the
        // energy method already sums shrinks as ctt/(N-1), while this term GROWS with
        // turn count (each added turn couples to the same core), which is what measured
        // toroid self-resonances demand. Returns ENERGY (J at the given per-turn
        // potentials), to be added to the self-pair energy before the 2E/dV^2 reduction.
        static double calculate_winding_to_core_self_energy(Coil coil, Core core,
                                                            const std::string& windingName,
                                                            const std::vector<double>& voltagesPerTurn,
                                                            std::optional<double> frequency = std::nullopt);

        // ABT #848: how much of an image plane the core is for the turns, from its MAS
        // permittivity (complex, with conduction) against the dielectric on its surface:
        // beta = (|eps_core| - eps_ext) / (|eps_core| + eps_ext). 1 for MnZn, nanocrystalline
        // and any conductor; ~0.6 for NiZn at 10 MHz; 1 when the database has neither
        // permittivity nor resistivity for the material. The floating-core terms above are
        // scaled by it when a frequency is given; with no frequency they keep beta = 1.
        static double core_image_factor(const Core& core, double frequency);

        std::map<std::pair<size_t, size_t>, double> calculate_capacitance_among_turns(Coil coil);

        // The optional core supplies the through-core inter-winding capacitance for
        // separated (non-adjacent) windings; omit it to keep the legacy behaviour where
        // separated windings have zero mutual capacitance.
        // frequency: where the capacitance is wanted (the impedance path passes its resonance);
        // it only sets the core image factor (see core_image_factor) — omit it for beta = 1.
        StrayCapacitanceOutput calculate_capacitance(Coil coil, std::optional<Core> core = std::nullopt, std::optional<double> frequency = std::nullopt);
        StrayCapacitanceOutput calculate_capacitance(Coil coil, OperatingPoint operatingPoint, std::optional<Core> core = std::nullopt, std::optional<double> frequency = std::nullopt);
    
    // Bipolar coordinate system for round-round energy density computation
    struct BipolarParams {
        double focalHalfDistance;
        double tau1;
        double tau2;
        double cosAngle;
        double sinAngle;
        double midX, midY;
    };
    static BipolarParams compute_bipolar_params(const Turn& t1, const Turn& t2);
    static double bipolar_tau_at_point(double lx, double ly, double a);
    static double bipolar_energy_density_at_point(double px, double py, const BipolarParams& params, double voltageDrop, double epsilonEff);

    // Internal three-input-multipole matrix used for the floating-node (V3) convergence and exposed
    // as the informational per-pair capacitance matrix. NOT the terminal Maxwell matrix (that is
    // calculate_maxwell_capacitance_matrix, built from the positive static capacitances).
    static ScalarMatrixAtFrequency calculate_capacitance_matrix_between_windings(double energy, double voltageDrop, double relativeTurnsRatio);
    // Canonical Biela/Kolar six-capacitor two-port network for a winding pair, built from the
    // static inter-winding capacitance C0 (energy method). Contains the intrinsic negative self
    // terms (-C0/6) — a valid math equivalent, NOT a passive netlist (see the .cpp for the node
    // topology). For circuit simulation use the positive tripole/pi-model instead.
    static SixCapacitorNetworkPerWinding calculate_six_capacitor_network(double staticInterwindingCapacitance);
        static std::vector<ScalarMatrixAtFrequency> calculate_maxwell_capacitance_matrix(Coil coil, std::map<std::string, std::map<std::string, double>> capacitanceAmongWindings);

};



class StrayCapacitanceOneLayer{

    public:

        StrayCapacitanceOneLayer(){
        };
        virtual ~StrayCapacitanceOneLayer() = default;
        // The core is optional for backwards compatibility, but for a TOROID it changes the
        // physics: the winding wraps the core cross-section (turn length from the core's
        // column width and depth, not a circle of the radial half-thickness), and the core is
        // a FLOATING conductor, so the turn-to-core capacitances form an energy-weighted
        // network that grows with the number of turns instead of the grounded-shield
        // cas/cab ladder, which converges to a turn-count-independent fixed point (ABT #845).
        double calculate_capacitance(Coil coil, std::optional<Core> core = std::nullopt);


};

// Effective relative permittivity of a wire's insulation/coating, resolved by wire type
// (ROUND: Albach enamel empirical formula; LITZ: serving/strand; FOIL/RECTANGULAR:
// coating material with a typical-film fallback; PLANAR: typical FR4). Exposed so FEM
// consumers (e.g. OMFEM) that mesh the coating as a dielectric can source the same value
// MKF uses, instead of re-deriving it. Definition in StrayCapacitance.cpp.
double get_wire_insulation_relative_permittivity(Wire wire);

} // namespace OpenMagnetics
