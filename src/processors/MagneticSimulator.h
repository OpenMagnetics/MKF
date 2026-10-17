#pragma once
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/CoreLosses.h"
#include "physical_models/WindingLosses.h"
#include "support/Utils.h"
#include "support/Settings.h"
#include "constructive_models/Mas.h"
#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {

// ABT #838: which winding-loss model the loss <-> temperature solve evaluates for a wound coil.
// The advisers' loss filters each publish their own winding-loss model (MagneticFilterLosses the
// full one, MagneticFilterLossesNoProximity DC + skin effect, the core adviser's post-processing
// DC only); the solve heats the network with exactly the winding losses the caller publishes, so
// the core and winding losses it returns are one self-consistent pair and are computed once.
enum class WindingLossesScope {
    FULL,                       // WindingLosses: DC + skin + proximity (simulate())
    OHMIC_AND_SKIN_EFFECT,      // DC + skin effect, no proximity
    OHMIC                       // DC only
};

// ABT #838: within one adviser run the same candidate is often solved again at the same
// operating point (by a filter and again by a later stage or flow). While a LossSolveMemoScope is
// alive on a thread, MagneticSimulator::calculate_core_losses keeps each converged solve of that
// thread under its exact inputs -- the serialized magnetic and operating point, the winding-loss
// scope, the simulator's models and every Settings value -- and answers an exact repeat with the
// stored result, which is what the solve returns. Nothing is matched approximately; a solve that
// throws is not stored. Scopes nest: an inner scope uses the outer scope's memo. The advisers open
// one per run; the memo is freed when the outermost scope ends.
// One converged loss <-> temperature solve: what calculate_core_losses returns and writes.
struct SolvedLosses {
    CoreLossesOutput coreLosses;
    OperatingPointExcitation excitation;
    std::optional<WindingLossesOutput> windingLosses;
    TemperatureOutput temperature;
    // The DC point of the magnetizing model at the solved core temperature (MagnetizingInductance::
    // get_last_dc_bias_point), for the saturation check's second criterion (ABT #1855). Empty when the
    // model ran no DC-bias fixed point, or magnetized a core held at its Curie temperature as paramagnetic.
    std::optional<DcBiasPoint> dcBiasPoint;
};

// ABT #1876 (Alf, 2026-10-14): the losses of an operating point evaluated at its ambient temperature,
// the core and the winding both, as the engine did before ABT #838: the loss <-> temperature solve's
// first iterate, before any temperature feedback.
struct LossesAtAmbient {
    CoreLossesOutput coreLosses;
    OperatingPointExcitation excitation;
    std::optional<WindingLossesOutput> windingLosses;
};

// ABT #1876 (Alf, 2026-10-14): the loss <-> temperature solve found that this operating point has NO
// thermal operating point. Exactly two cases:
//   - the thermal network did not converge on the losses of a loop iterate
//     (ThermalNetworkNotConvergedException; code CALCULATION_INVALID_RESULT), or
//   - the winding runs away with the core held at its data ceiling (code CALCULATION_DIVERGED).
// The message is MKF's, with the loop's history. It carries the point's losses at ambient, which
// MagneticSimulator::simulate publishes for that point (and no temperature); every other caller sees it
// as the CalculationException it is. Any other error of the solve is not this exception.
class NoThermalOperatingPointException : public CalculationException {
  public:
    NoThermalOperatingPointException(ErrorCode code, const std::string& message, LossesAtAmbient lossesAtAmbient)
        : CalculationException(code, message), _lossesAtAmbient(std::move(lossesAtAmbient)) {}
    const LossesAtAmbient& get_losses_at_ambient() const { return _lossesAtAmbient; }
  private:
    LossesAtAmbient _lossesAtAmbient;
};

// ABT #1876 (Alf, 2026-10-14): what MagneticSimulator::calculate_core_losses does for an operating point with no
// thermal operating point (NoThermalOperatingPointException).
enum class NoThermalOperatingPointPolicy {
    // The public per-operating-point answer, as simulate gives it per point: the core (and, when asked, winding)
    // losses evaluated at ambient (LossesAtAmbient) and NO temperature output (*temperatureOutput left empty).
    // The returned CoreLossesOutput.temperature is that ambient, the temperature the losses were evaluated at.
    LOSSES_AT_AMBIENT,
    // The advisers and ResultsAnalyzer: the NoThermalOperatingPointException propagates, as it always did.
    THROW
};

// ABT #838 / #1855 / #1908: the highest temperature up to which every temperature-dependent material
// datum the loss <-> temperature solve evaluates the core with still answers, and which datum ends there.
// The solve keeps the core below it, and holds it there when the losses drive it to it.
enum class CoreDataCeilingSource {
    CURIE_TEMPERATURE,                  // the material's Curie temperature: above it the core is paramagnetic
    PERMEABILITY_TABLE_END,             // the last point of a mu_i(T) table that ends below Curie (ABT #1739)
    STEINMETZ_TEMPERATURE_FACTOR_END    // where an evaluated Steinmetz range's ct(T) stops being positive
};
struct CoreDataCeiling {
    double temperature;
    CoreDataCeilingSource source;
};

class LossSolveMemoScope {
    public:
        LossSolveMemoScope();
        ~LossSolveMemoScope();
        LossSolveMemoScope(const LossSolveMemoScope&) = delete;
        LossSolveMemoScope& operator=(const LossSolveMemoScope&) = delete;
        struct Statistics {
            size_t hits = 0;
            size_t misses = 0;
            size_t entries = 0;
        };
        // The memo active on this thread (all zero when no scope is open).
        static Statistics get_statistics();
    private:
        bool _owner = false;
};

class MagneticSimulator {
    private:
        std::optional<DcBiasPoint> _lastDcBiasPoint;
        ReluctanceModels _reluctanceModelName;

        MagnetizingInductance _magnetizingInductanceModel;
        CoreLosses _coreLossesModel;

        SolvedLosses solve_losses_at_operating_temperature(const OperatingPoint& operatingPoint, Magnetic& magnetic,
                                                           WindingLossesScope windingLossesScope);
        std::string loss_solve_key(const OperatingPoint& operatingPoint, const Magnetic& magnetic,
                                   WindingLossesScope windingLossesScope) const;

    public:
        // The magnetic with its core material made paramagnetic (mu_r = 1, its own material name): what the
        // magnetizing model sees of a core the loss solve holds at its Curie temperature (ABT #1855).
        static Magnetic make_paramagnetic_magnetic(const Magnetic& magnetic);

        MagneticSimulator() {
            // Use Settings instead of defaults
            auto& settings = Settings::GetInstance();
            _reluctanceModelName = settings.get_reluctance_model();

            _magnetizingInductanceModel = MagnetizingInductance(to_string(_reluctanceModelName));
        }

        void set_reluctance_model_name(ReluctanceModels model) {
            _reluctanceModelName = model;
        }
        // Default preference: heads the core-loss model cascade (see CoreLosses::set_core_losses_model_name).
        void set_core_losses_model_name(CoreLossesModels model) {
            _coreLossesModel.set_core_losses_model_name(model);
        }
        // Explicit request: exactly this model, throwing for a material that cannot run it (ABT #1497).
        void set_core_losses_requested_model_name(CoreLossesModels model) {
            _coreLossesModel.set_core_losses_requested_model_name(model);
        }

        // The design's topology for the magnetizing-inductance model (see
        // MagnetizingInductance::set_topology). simulate() sets it from its Inputs; callers of the
        // per-operating-point methods below set it from theirs.
        void set_topology(std::optional<Topology> topology) {
            _magnetizingInductanceModel.set_topology(topology);
        }

        // One output per operating point. A point with no thermal operating point
        // (NoThermalOperatingPointException, ABT #1876, Alf 2026-10-14) gets its core and winding losses
        // evaluated at ambient (LossesAtAmbient) and NO temperature output; the other points are simulated
        // normally, and ResultsAnalyzer::analyze_results reports NO_THERMAL_OPERATING_POINT for it. Every
        // other error still throws for the whole design. fastMode publishes no temperature for ANY point
        // (and skips the leakage inductance); analyze_results tells the two apart by solving again.
        Mas simulate(Mas mas, bool fastMode=false);
        Mas simulate(const Inputs& inputs, const Magnetic& magnetic, bool fastMode=false);

        // Builds a manufacturer-style datasheet (MagneticManufacturerInfo.datasheetInfo) from an
        // already-simulated MAS, reading the per-operating-point Outputs and the magnetic geometry.
        // Existing manufacturerInfo fields (name, reference, description, family, ...) are preserved;
        // only the datasheetInfo block is populated/overwritten. The result is also attached to
        // mas.magnetic.manufacturerInfo. Requires mas to carry simulation outputs (run simulate first).
        MagneticManufacturerInfo build_datasheet(Mas& mas);
        // Full thermal-network hot-spot for one simulated operating point, configured through
        // TemperatureConfig::fromSimulatedOutput (ABT #906) — the same config the temperature-map
        // plot wrappers use, so the exported MAS and the UI can never disagree.
        TemperatureOutput calculate_temperature(OperatingPoint& operatingPoint, Magnetic magnetic, const Outputs& output);
        // ABT #838: core losses at the temperature the part's dissipation produces, by
        // fixed-point iteration from ambient on the thermal network (the core-loss <-> temperature
        // loop this function has always held, now always run to convergence):
        //   - coil WOUND (turns description present): core AND winding heat the same full network,
        //     so neither temperature can be found without the other's losses. Core losses are
        //     evaluated at T_core = hottest CORE node, winding losses at T_winding = mean TURN-node
        //     temperature (WindingLosses takes one temperature for the whole coil). simulate() uses
        //     this, and so do the advisers' loss filters and the PyOM/WASM wrappers when they pass a
        //     wound coil.
        //   - no turns (the core advisers' sentinel-coil candidates): no winding losses exist to
        //     compute; the core-only network carries the core losses alone.
        // Converged when the network returns both temperatures within
        // kOperatingTemperatureTolerance of the ones the losses were evaluated at. Each
        // temperature's step is halved whenever its residual changes sign (losses falling steeply
        // with temperature make plain substitution oscillate), and the core is kept below the
        // highest temperature where its material is ferromagnetic and its mu_i(T) has data
        // (InitialPermeability::get_maximum_ferromagnetic_data_temperature: the last table point
        // when the table ends below Curie, else Curie), and where the evaluated Steinmetz ranges' ct(T)
        // are positive (get_core_data_ceiling, ABT #1908). Not settling within
        // kOperatingTemperatureMaximumIterations THROWS (there is no operating point to report, and
        // the last iterate is not one). Being driven to that ceiling holds the core there and settles the
        // winding alone; the held core is magnetized with mu_r = 1 only when the ceiling is its Curie
        // temperature, and keeps its last tabulated permeability when it is the end of a mu_i(T) table
        // below Curie (explicit exception approved by Alf, 2026-10-14,
        // ABT #1855; analyze_results reports it). The thermal network failing to converge on an
        // iterate's losses, and the held winding running away, throw NoThermalOperatingPointException
        // carrying the losses at ambient (ABT #1876).
        // `windingLossesOutput` / `temperatureOutput`, when given, receive the winding losses and
        // the full-network temperature solved together with the returned core losses (asking for
        // winding losses of an unwound coil throws). `windingLossesScope` selects the winding-loss
        // model the solve evaluates (and returns) for a wound coil.
        // `noThermalOperatingPointPolicy` (ABT #1876, Alf 2026-10-14): by default a point with no thermal
        // operating point returns its losses at ambient and leaves `temperatureOutput` EMPTY (the public
        // per-point entry point, as simulate does per point); THROW (the advisers, ResultsAnalyzer) lets the
        // NoThermalOperatingPointException propagate.
        // K, both temperatures. 2 K (was 0.5 K; approved by Alf 2026-10-11): on the 522 solves of the
        // web9 magnetic adviser run it moves the hot spot by <= 0.57 K, the winding temperature by
        // <= 2 K and the total loss by <= 0.96%, and saves ~300 of 1775 iterations.
        static constexpr double kOperatingTemperatureTolerance = 2.0;
        static constexpr size_t kOperatingTemperatureMaximumIterations = 100;
        // ABT #1855: whether a solve held the core at its material's data ceiling (the losses drove it
        // there, solve_losses_at_operating_temperature): its core losses were evaluated within
        // kOperatingTemperatureTolerance below InitialPermeability::get_maximum_ferromagnetic_data_temperature
        // and the network's maximum temperature lies above that ceiling. This overload knows only that
        // mu_i(T) / Curie part (get_ferromagnetic_data_ceiling); the full ceiling the solve uses also depends
        // on the excitation (get_core_data_ceiling) and is judged by the overload below.
        static bool is_core_held_at_data_ceiling(const CoreMaterial& coreMaterial, double coreLossesTemperature,
                                                 double maximumTemperature);
        // The same against a ceiling of get_core_data_ceiling (none: never held).
        static bool is_core_held_at_data_ceiling(const std::optional<CoreDataCeiling>& ceiling, double coreLossesTemperature,
                                                 double maximumTemperature);
        // The mu_i(T) / Curie part of the ceiling: InitialPermeability::get_maximum_ferromagnetic_data_temperature,
        // CURIE_TEMPERATURE when it is the Curie temperature, else PERMEABILITY_TABLE_END. None without either.
        static std::optional<CoreDataCeiling> get_ferromagnetic_data_ceiling(const CoreMaterial& coreMaterial);
        // ABT #1908: the ceiling the loss <-> temperature solve keeps the core below, for `core` priced by
        // `coreLossesModel` on `excitation` (the excitation the solve passes to calculate_core_losses) from
        // `ambientTemperature` up: the lower of get_ferromagnetic_data_ceiling and the end of positive ct(T)
        // of the Steinmetz ranges the losses evaluate (CoreLosses::get_temperature_factor_data_end). The
        // saturation field strength sets no ceiling of its own: past its data it goes linearly to zero at the
        // Curie temperature (Core::get_magnetic_field_strength_saturation), which the first part already is. On a tie the earlier source in
        // CoreDataCeilingSource names it. Every part is derived from the material's data.
        static std::optional<CoreDataCeiling> get_core_data_ceiling(const Core& core, const OperatingPointExcitation& excitation,
                                                                    double ambientTemperature, CoreLosses& coreLossesModel);
        // The DC point of the last calculate_core_losses' solve (SolvedLosses::dcBiasPoint).
        std::optional<DcBiasPoint> get_last_dc_bias_point() const { return _lastDcBiasPoint; }
        CoreLossesOutput calculate_core_losses(OperatingPoint& operatingPoint, Magnetic magnetic,
                                               std::optional<WindingLossesOutput>* windingLossesOutput = nullptr,
                                               std::optional<TemperatureOutput>* temperatureOutput = nullptr,
                                               WindingLossesScope windingLossesScope = WindingLossesScope::FULL,
                                               NoThermalOperatingPointPolicy noThermalOperatingPointPolicy =
                                                   NoThermalOperatingPointPolicy::LOSSES_AT_AMBIENT);
        LeakageInductanceOutput calculate_leakage_inductance(OperatingPoint& operatingPoint, Magnetic magnetic);
        static LeakageInductanceOutput calculate_leakage_inductance(Magnetic magnetic, double frequency);
        MagnetizingInductanceOutput calculate_magnetizing_inductance(OperatingPoint& operatingPoint, Magnetic magnetic);
        WindingLossesOutput calculate_winding_losses(OperatingPoint& operatingPoint, Magnetic magnetic, std::optional<double> temperature = std::nullopt);

};


} // namespace OpenMagnetics