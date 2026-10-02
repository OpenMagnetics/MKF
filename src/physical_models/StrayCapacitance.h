#pragma once
#include <array>
#include <complex>
#include <functional>
#include <set>
#include <limits>
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
        // insulationLayersInBetween: when given, used instead of looking the layers up from the turns'
        // coordinates (the toroid pair split evaluates the outer crossings, whose coordinates are
        // not in the frame the layer lookup works in; the layers between two turns are the same
        // on every run of the ring).
        std::vector<double> preprocess_data_for_round_wires(Turn firstTurn, Wire firstWire, Turn secondTurn, Wire secondWire, std::optional<Coil> coil = std::nullopt, std::optional<std::vector<Layer>> insulationLayersInBetween = std::nullopt);
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
        // WHICH MODEL ACTUALLY COMPUTED THE PAIRS, which is not always the one selected. Any pair
        // involving a flat conductor is routed to ParallelPlate regardless of the caller's choice,
        // so reporting the SELECTED name for such a winding says the model ran when it did not.
        // Recorded per pair as the dispatch decides, and reported at the end (ABT #950).
        std::set<std::string> _methodsUsed;
        static double calculate_area_between_two_turns_using_diagonals(Turn firstTurn, Turn secondTurn);
        static double calculate_area_between_two_turns_using_vecticals_and_horizontals(Turn firstTurn, Turn secondTurn);
        StrayCapacitanceOutput calculate_capacitance_with_voltages(Coil coil, std::map<std::string, double> voltageRmsPerWinding, std::optional<Core> core = std::nullopt, std::optional<double> frequency = std::nullopt, std::optional<CoreElectricalReference> coreElectricalReference = std::nullopt);
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
        //
        // Neighbour screening: the element is the turn's PARTIAL capacitance to the core, i.e. the
        // row sum of the Maxwell capacitance matrix -- the charge on the turn when it and every
        // neighbouring conductor sit at the same potential. The field a turn sends towards the
        // core between itself and a neighbour ends on that neighbour instead, and that part is
        // the turn-to-turn element's business. leftNeighbourPitch / rightNeighbourPitch are the
        // centre-to-centre distances, along the core surface, to the nearest conductor on each
        // side that faces the same surface; infinity (the default) means no neighbour on that
        // side and reproduces the isolated cylinder over a plane. Each side screens half the
        // turn: element = (C'(p_left) + C'(p_right)) / 2 * length, with C'(p) the per-length
        // capacitance of one conductor of an infinite periodic row of pitch p over the plane
        // (calculate_conductor_row_over_plane_capacitance_per_length).
        static double calculate_turn_to_core_capacitance(double conductingRadius, double turnLength,
                                                         double wireCoatingThickness, double wireCoatingRelativePermittivity,
                                                         double airGapToCore,
                                                         double coreCoatingThickness, double coreCoatingRelativePermittivity,
                                                         double bobbinThickness = 0.0, double bobbinRelativePermittivity = 1.0,
                                                         double leftNeighbourPitch = std::numeric_limits<double>::infinity(),
                                                         double rightNeighbourPitch = std::numeric_limits<double>::infinity());

        // Per-unit-length capacitance to a grounded conducting plane of ONE conductor of an
        // infinite periodic row of equipotential circular conductors (radius r, axis at height
        // H > r above the plane, centre-to-centre pitch p >= 2r along the plane). Solved to
        // ~1e-5 by the charge simulation method with the exact periodic Green's function of
        // the row and its image. Limits: p -> infinity gives the isolated cylinder over a plane,
        // 2 pi eps0 / acosh(H/r) (Smythe), returned in closed form for p = infinity; r << p, H
        // gives the thin-wire grid 2 pi eps0 / ln(sinh(2 pi H/p) / sinh(pi r/p)). Throws on
        // H <= r, p < 2r, or a solution that does not meet the boundary condition.
        static double calculate_conductor_row_over_plane_capacitance_per_length(double radius, double axisHeight, double pitch);

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

        // ABT #1165: the energy stored in the turn-to-core elements of a winding PAIR against
        // the floating core -- the turn -> core -> turn path, which exists for every pair and
        // not only for separated ones. The first winding's turns sit at
        // firstWindingPotentialOffset + V_i (the offset is the orchestration loop's V3, the
        // common-mode offset it already applies to the turn-to-turn drops) and the second
        // winding's at -V_j (opposing DM currents, the sign flip the loop already applies);
        // the core is one floating node whose potential balances the charge over BOTH windings.
        // Returns ENERGY (J) so the caller folds it into the same sum as the turn-to-turn pairs
        // and reduces once, instead of reducing a second capacitance against a second voltage.
        static double calculate_winding_pair_to_core_energy(Coil coil, Core core,
                                                            const std::string& firstWindingName,
                                                            const std::string& secondWindingName,
                                                            const std::vector<double>& voltagesPerTurn,
                                                            double firstWindingPotentialOffset = 0.0,
                                                            std::optional<double> frequency = std::nullopt);

        // ABT #1167: the potential the core node is held at, or nullopt when it FLOATS (its
        // potential set by charge balance over the turns facing it, which is what every caller
        // assumed before magnetic.coreElectricalReference existed). An ABSENT reference, and an
        // explicit "floating" one, both return nullopt -- bit-for-bit today's behaviour.
        //  - grounded: a node with no potential swing, i.e. 0 in the frame the per-turn
        //    potentials are expressed in (every winding's "end" terminal sits at 0 there, so
        //    every isolation side's local ground is the same 0; isolationSide names the node,
        //    it does not move it).
        //  - tiedToWinding: the potential of the named winding's start or end TERMINAL, taken
        //    from the voltage dividers (exact terminal potentials, where voltagesPerTurn carries
        //    turn CENTRES). Throws if the named winding is not in the coil, or if the required
        //    winding/terminal fields are missing -- never a silent fall back to floating.
        static std::optional<double> resolve_core_reference_potential(
                Coil& coil,
                const std::optional<CoreElectricalReference>& coreElectricalReference,
                const StrayCapacitanceOutput& voltagesOutput,
                const std::map<std::string, double>& voltageRmsPerWinding);

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
                                                            std::optional<double> frequency = std::nullopt,
                                                            // ABT #1167: nullopt = FLOATING core (charge-balanced node, the
                                                            // pre-existing behaviour); a value pins the core at that potential,
                                                            // which is what a clip, strap or flux band to a circuit node does.
                                                            std::optional<double> fixedCorePotential = std::nullopt);

        // ABT #848: how much of an image plane the core is for the turns, from its MAS
        // permittivity (complex, with conduction) against the dielectric on its surface:
        // beta = (|eps_core| - eps_ext) / (|eps_core| + eps_ext). 1 for MnZn, nanocrystalline
        // and any conductor; ~0.6 for NiZn at 10 MHz; 1 when the database has neither
        // permittivity nor resistivity for the material. The floating-core terms above are
        // scaled by it when a frequency is given; with no frequency they keep beta = 1.
        static double core_image_factor(const Core& core, double frequency);

        // ===================== ABT #1166 BEGIN (air gap in the through-core path) =====================
        // The through-core path (turn -> core -> turn) above treats the ferrite as ONE floating
        // equipotential node. An air gap breaks metal-to-metal continuity and inserts
        // eps0*eps_r*A/g into that path -- but ONLY where the gap cuts EVERY conductive route
        // between the two windings' footprints. Design note "Transformer Stray Capacitance"
        // (2026-09-09) §13-§15 classifies the three cases, all decidable from MAS data:
        //
        //   A  centre-leg gap only .................. the outer legs still touch, the core is still
        //                                             one body                     -> SHARED_CORE_NODE
        //   B  all legs gapped, CONCENTRIC windings .. each winding faces BOTH halves and those
        //                                             couplings shunt the gap      -> SHARED_CORE_NODE
        //   C  all legs gapped AND the two windings sit on OPPOSITE sides of the gap plane
        //      (side-by-side / split bobbin, each facing one half) .............. -> SPLIT_CORE_NODES
        //
        // In case C the gap capacitance is IN SERIES with Cpc and Csc and, being the smallest
        // element, sets the total: the common offline-flyback split-bobbin construction, whose
        // inter-winding term drives its common-mode noise model, so the single-node answer is
        // about an order of magnitude high at a 1 mm gap.
        //
        // What the gap does NOT change: the INTRA-winding (self) term. Each core half is still
        // locally equipotential below the ferrite's dielectric relaxation frequency, so a
        // winding's own terminal shunt through the core -- calculate_winding_to_core_self_energy
        // -- is unaffected, and it is deliberately left alone. The gap moves the common-mode path
        // BETWEEN the windings only.
        //
        // Powder / distributed-gap materials (MPP, Kool Mu, iron powder), like NiZn, are
        // high-resistivity THROUGHOUT, so the single-node picture fails for them for a different
        // reason -- imaging, not continuity. That is core_image_factor's job and it is not
        // re-done here: a distributed-gap core carries no discrete non-residual gapping, so it
        // never reaches SPLIT_CORE_NODES and the two effects cannot double-count.
        enum class ThroughCoreGapTopology { SHARED_CORE_NODE, SPLIT_CORE_NODES };

        struct ThroughCoreGapSplit {
            ThroughCoreGapTopology topology = ThroughCoreGapTopology::SHARED_CORE_NODE;
            // Series capacitance joining the two core bodies, F. Only meaningful (and only
            // strictly positive) when topology == SPLIT_CORE_NODES.
            double gapCapacitance = 0;
            // Axial coordinate of the common gap plane, in the core-centred frame the turn and
            // bobbin coordinates use. Only meaningful when topology == SPLIT_CORE_NODES.
            double gapPlaneAxialCoordinate = 0;
            // Total gapped cross-section summed over the columns, m^2, and the (single) gap
            // length, m -- reported for diagnostics and for the tests.
            double totalGappedArea = 0;
            double gapLength = 0;
        };

        // Cgap = fringingFactor * eps0 * eps_r * A / g, with A the TOTAL gapped cross-section
        // across all the gapped columns (the columns' gaps are in parallel between the two core
        // bodies) and g the gap length. THROWS on a non-positive area, length or permittivity --
        // a gap whose area or length the record does not carry is a missing input, not a number
        // to invent.
        static double gap_capacitance(double totalGappedArea, double gapLength,
                                      double gapRelativePermittivity, double fringingFactor);

        // Classifies the core+coil into case A/B (SHARED_CORE_NODE, the historical behaviour) or
        // case C (SPLIT_CORE_NODES, with the gap capacitance filled in). Purely geometric and
        // free of the capacitance model, so it can be exercised on its own.
        static ThroughCoreGapSplit core_gap_topology(Core core, Coil coil,
                                                     const std::string& firstWindingName,
                                                     const std::string& secondWindingName);
        // The same classification for a caller that already holds the core, the wound coil and its
        // turns (the through-core path): nothing is copied again. windingTurns may be null, in which
        // case the coil is wound if needed and its turns read, exactly as the by-value form does.
        // Identical result either way; this exists because the classification runs inside every
        // through-core evaluation, and copying the coil and its turns twice more per call cost up
        // to 57% of that call's runtime (ABT #1200).
        static ThroughCoreGapSplit core_gap_topology(Core& core, Coil& coil,
                                                     const std::vector<Turn>* windingTurns,
                                                     const std::string& firstWindingName,
                                                     const std::string& secondWindingName);
        // ====================== ABT #1166 END ======================

        std::map<std::pair<size_t, size_t>, double> calculate_capacitance_among_turns(Coil coil);

        // Capacitance seen at the DIFFERENTIAL-MODE port of a two-winding magnetic: the two
        // windings in series opposition, far ends joined, driven across their two start
        // terminals. Energy method over EVERY capacitance the port charges -- the turn-to-turn
        // pairs inside each winding and between them, and the turn-to-core elements against the
        // floating (or bonded) core -- with the DM potentials (first winding +v_i, second -v_j,
        // joined far ends at 0): C_DM = 2 W / V_port^2. This is what shunts the leakage
        // inductance at the DM port; the inter-winding entry of the capacitance matrix is one
        // branch of it only.
        double calculate_differential_mode_capacitance(Coil coil, Core core, std::optional<double> frequency = std::nullopt,
                                                       std::optional<CoreElectricalReference> coreElectricalReference = std::nullopt);

        // The optional core supplies the through-core inter-winding capacitance for
        // separated (non-adjacent) windings; omit it to keep the legacy behaviour where
        // separated windings have zero mutual capacitance.
        // frequency: where the capacitance is wanted (the impedance path passes its resonance);
        // it only sets the core image factor (see core_image_factor) — omit it for beta = 1.
        StrayCapacitanceOutput calculate_capacitance(Coil coil, std::optional<Core> core = std::nullopt, std::optional<double> frequency = std::nullopt, std::optional<CoreElectricalReference> coreElectricalReference = std::nullopt);
        StrayCapacitanceOutput calculate_capacitance(Coil coil, OperatingPoint operatingPoint, std::optional<Core> core = std::nullopt, std::optional<double> frequency = std::nullopt, std::optional<CoreElectricalReference> coreElectricalReference = std::nullopt);
        // The core's electrical reference is a property of the MAGNETIC, not of the core or the
        // coil (it describes how the assembled part is bonded), so these overloads are the ones
        // that can honour it without the caller restating it: they read
        // magnetic.coreElectricalReference. The Coil+Core overloads above stay for callers that
        // hold no magnetic; an absent reference there means floating, which is exactly what they
        // computed before ABT #1167.
        StrayCapacitanceOutput calculate_capacitance(const Magnetic& magnetic, std::optional<double> frequency = std::nullopt);
        StrayCapacitanceOutput calculate_capacitance(const Magnetic& magnetic, OperatingPoint operatingPoint, std::optional<double> frequency = std::nullopt);
    
    // ---- The model's electrostatic energy, element by element (for the Painter) -------------------------
    // Every capacitance this class reports is 2W/dV^2 of an energy W that is a sum over elements: the
    // turn-to-turn pairs (Albach/Koch/... for round conductors, evaluated per crossing on a toroid; parallel
    // plates for flat ones) and the screened turn-to-core faces against the core node.
    // calculate_electric_energy_elements lists those elements, built by the same code paths as
    // calculate_capacitance_among_turns and turn_to_core_elements (same pairs, same toroid split, same faces,
    // same screening, same image factor), with the geometry each one is evaluated at, so that a field picture
    // can put each element's energy where the element's own field puts it.
    struct ElectricEnergyElement {
        enum class Kind { ROUND_PAIR, PLATE_PAIR, CORE_FACE };
        // The core surface a CORE_FACE element faces. NONE: the face builder had no bobbin geometry to place
        // the turn against; such an element has a capacitance but no position, and sampling it throws.
        enum class Surface { NONE, BORE, OUTSIDE, TOROID_FLAT_FACE, COLUMN, UPPER_YOKE, LOWER_YOKE };
        Kind kind;
        size_t firstTurnIndex = 0;
        size_t secondTurnIndex = 0;                 // pairs only
        std::array<double, 2> firstCentre = {0, 0}; // the crossing (or turn) centre the element is evaluated at, m
        std::array<double, 2> secondCentre = {0, 0};
        double length = 0;                          // m, the length share this element stands for
        double capacitance = 0;                     // F, the element over its length (core image factor included)
        bool inPlane = true;                        // false for a toroid's top and bottom runs (not in the ring plane)
        // ROUND_PAIR (Albach's flux tubes): conductor radius r0, coating thickness delta and permittivity, the
        // gap between the coating surfaces and its effective permittivity (as the pair model computed them)
        double conductingRadius = 0;
        double coatingThickness = 0;
        double coatingPermittivity = 1;
        double gap = 0;
        double gapPermittivity = 1;
        // PLATE_PAIR: the gap rectangle between the facing conductor surfaces
        std::array<double, 2> plateCentre = {0, 0};
        std::array<double, 2> plateHalfExtents = {0, 0};
        // CORE_FACE: conductingRadius is the equivalent radius of the turn-to-core element; airEquivalentGap is
        // sum t_i/eps_i of the stack (coating, air, core jacket, bobbin wall); realDistance is the distance from
        // the conductor axis to the ferrite surface in the real geometry (outer half extent + air + bobbin wall
        // + core jacket); the pitches are the row pitches of the element, positivePitch on the side the surface
        // tangent points to (rectangular column: +y; yokes: +x; toroid: counter-clockwise).
        Surface surface = Surface::NONE;
        double airEquivalentGap = 0;
        double realDistance = 0;
        double positivePitch = std::numeric_limits<double>::infinity();
        double negativePitch = std::numeric_limits<double>::infinity();
    };
  private:
    // Set only while calculate_electric_energy_elements runs calculate_capacitance_among_turns: every pair
    // element the turn-to-turn dispatch evaluates is appended here with its geometry.
    std::vector<ElectricEnergyElement>* _energyElementRecorder = nullptr;
  public:
    std::vector<ElectricEnergyElement> calculate_electric_energy_elements(Coil coil, std::optional<Core> core = std::nullopt, std::optional<double> frequency = std::nullopt);

    // The energy each element stores (J) at the given complex rms turn potentials (indexed like the turns
    // description), 1/2 C |dV|^2. The core is held at fixedCorePotential when given (a bonded core); otherwise it
    // floats: one charge-balanced node, sum C_i (V_i - V_c) = 0 over all its faces, or, for two windings on a core
    // that core_gap_topology classifies SPLIT_CORE_NODES, two bodies joined by the gap capacitance (each winding
    // on its own body), whose energy is returned separately as gapEnergy.
    struct ElectricEnergyDistribution {
        std::vector<double> energyPerElement;
        std::vector<std::complex<double>> corePotentialPerBody;
        double gapEnergy = 0;
    };
    static ElectricEnergyDistribution calculate_electric_energy_per_element(const std::vector<ElectricEnergyElement>& elements, Coil coil, std::optional<Core> core,
                                                                             const std::vector<std::complex<double>>& turnPotentials,
                                                                             std::optional<std::complex<double>> fixedCorePotential = std::nullopt);

    // The charge-simulation solution behind calculate_conductor_row_over_plane_capacitance_per_length, in
    // lengths normalised to the radius and unit conductor potential: phi(x, y) = sum q_k G(x, y; x_k, y_k) with
    // the periodic row-and-image Green's function G. capacitancePerLength = 2 pi eps0 sum q_k.
    struct ConductorRowOverPlaneSolution {
        double height = 0;   // H/r
        double period = 0;   // p/r
        std::vector<double> sourceX;
        std::vector<double> sourceY;
        std::vector<double> charges;
        double capacitancePerLength = 0;
    };
    static const ConductorRowOverPlaneSolution& solve_conductor_row_over_plane(double radius, double axisHeight, double pitch);
    // grad phi (normalised) of a solution at (x, y)
    static std::array<double, 2> conductor_row_over_plane_potential_gradient(const ConductorRowOverPlaneSolution& solution, double x, double y);

    // Spatial distribution of one element's energy per unit length (J/m) in the plane, as weighted points
    // (x, y, J/m) at a spacing finer than samplingSpacing. Each element is sampled with its own field:
    //  - ROUND_PAIR: Albach's flux tubes (straight lines between the coating surfaces, radial in the coating,
    //    the coating taken to second order in delta/(r0 + delta) as in Albach's closed form), each tube a series
    //    capacitor eps0 R cos(phi) dphi / (2 L R cos(phi) + s_eq(phi)) carrying the pair's dV.
    //  - PLATE_PAIR: the uniform field of the parallel-plate element in its gap rectangle.
    //  - CORE_FACE: the field of the conductor (row) over the plane in the air-equivalent geometry, integrated
    //    in bipolar coordinates, with each half cell at its own pitch, carried to the real geometry (the
    //    air-equivalent height stretched over the real distance; wrapped around the ring for a toroid).
    // The points carry exactly energyPerLength between them. Returns the energy per length that has no place in
    // the plane (beyond half a toroid surface's circumference).
    using EnergySampleCallback = std::function<void(double x, double y, double energyPerLength)>;
    static double sample_electric_energy_element(const ElectricEnergyElement& element, double energyPerLength, double samplingSpacing,
                                                 const EnergySampleCallback& callback);
    // The energy per unit length, per V^2 of potential difference, that an element's own field integrates to
    // (the physics sample_electric_energy_element distributes): ROUND_PAIR the flux-tube integral, CORE_FACE the
    // bipolar integral of the (row) field. The closed-form element is capacitance / (2 length).
    static double integrate_electric_energy_element(const ElectricEnergyElement& element, double samplingSpacing);

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
