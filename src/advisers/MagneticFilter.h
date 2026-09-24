#pragma once
#include "physical_models/Temperature.h"
#include "physical_models/CoreLosses.h"
#include "physical_models/Impedance.h"
#include "physical_models/MagneticEnergy.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "constructive_models/Coil.h"
#include "processors/Inputs.h"
#include "processors/MagneticSimulator.h"
#include "processors/Outputs.h"
#include "constructive_models/Core.h"
#include "constructive_models/Mas.h"
#include "Definitions.h"
#include <cmath>
#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {

class MagneticFilter {
    public: 
        static std::shared_ptr<MagneticFilter> factory(MagneticFilters filterName, std::optional<Inputs> inputs = std::nullopt);

        MagneticFilter() { };
        virtual ~MagneticFilter() = default;
        virtual std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr) = 0;

        // Whether this filter can judge `magnetic` at all. Most filters compute from the core and
        // coil, so a datasheet-only catalogue part (no construction, see Magnetic.h) is outside
        // anything they can say: the catalogue adviser then records NO score for that part and
        // filter -- never a made-up one -- and ranks the part on the filters that do apply to it.
        // Filters that can read their answer from the datasheet override this.
        virtual bool applies_to(Magnetic* magnetic) const { return magnetic->has_core() && magnetic->has_coil(); }
};

class MagneticFilterAreaProduct : public MagneticFilter {
    private:
        std::map<std::string, double> _materialScaledMagneticFluxDensities;
        std::map<std::string, double> _bobbinFillingFactors;
        OperatingPointExcitation _operatingPointExcitation;
        std::vector<double> _areaProductRequiredPreCalculations;
        WindingSkinEffectLosses _windingSkinEffectLossesModel;
        std::shared_ptr<CoreLossesModel> _coreLossesModelSteinmetz;
        std::shared_ptr<CoreLossesModel> _coreLossesModelProprietary;
        double _averageMarginInWindingWindow = 0;
        double _magneticFluxDensityReference = 0.18;

    public:
        MagneticFilterAreaProduct() {};
        MagneticFilterAreaProduct(Inputs inputs);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        double get_estimated_area_product_required(Inputs inputs);
};

class MagneticFilterEnergyStored : public MagneticFilter {
    private:
        std::map<std::string, std::string> _models;
        MagneticEnergy _magneticEnergy;
        double _requiredMagneticEnergy;

    public:
        MagneticFilterEnergyStored() {};
        MagneticFilterEnergyStored(Inputs inputs);
        MagneticFilterEnergyStored(Inputs inputs, std::map<std::string, std::string> models);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterCost : public MagneticFilter {
    public:
        MagneticFilterCost() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterEstimatedCost : public MagneticFilter {
    private:
        double _estimatedParallels = 0;
        double _estimatedWireTotalArea = 0;
        double _wireAirFillingFactor = 0;
        double _skinDepth = 0;
        double _averageMarginInWindingWindow = 0;

    public:
        MagneticFilterEstimatedCost() {};
        MagneticFilterEstimatedCost(Inputs inputs);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterCoreAndDcLosses : public MagneticFilter {
    private:
        MagnetizingInductance _magnetizingInductance;
        WindingOhmicLosses _windingOhmicLosses;
        std::map<std::string, std::string> _models;
        std::shared_ptr<CoreLossesModel> _coreLossesModelSteinmetz = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Steinmetz"}}));
        std::shared_ptr<CoreLossesModel> _coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}));
        double _maximumPowerMean = 0;
    public:

        MagneticFilterCoreAndDcLosses();
        MagneticFilterCoreAndDcLosses(Inputs inputs);
        MagneticFilterCoreAndDcLosses(Inputs inputs, std::map<std::string, std::string> models);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterCoreDcAndSkinLosses : public MagneticFilter {
    private:
        MagnetizingInductance _magnetizingInductance;
        WindingOhmicLosses _windingOhmicLosses;
        WindingSkinEffectLosses _windingSkinEffectLosses;
        std::map<std::string, std::string> _models;
        std::shared_ptr<CoreLossesModel> _coreLossesModelSteinmetz = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Steinmetz"}}));
        std::shared_ptr<CoreLossesModel> _coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}));
        double _maximumPowerMean = 0;
    public:

        MagneticFilterCoreDcAndSkinLosses();
        MagneticFilterCoreDcAndSkinLosses(Inputs inputs);
        MagneticFilterCoreDcAndSkinLosses(Inputs inputs, std::map<std::string, std::string> models);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterLosses : public MagneticFilter {
    private:
        std::map<std::string, std::string> _models;
        MagneticSimulator _magneticSimulator;
    public:
        MagneticFilterLosses() {};
        MagneticFilterLosses(std::map<std::string, std::string> models);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterLossesNoProximity : public MagneticFilter {
    private:
        std::map<std::string, std::string> _models;
        WindingOhmicLosses _windingOhmicLosses;
        WindingSkinEffectLosses _windingSkinEffectLosses;
        MagneticSimulator _magneticSimulator;
    public:
        MagneticFilterLossesNoProximity() {};
        MagneticFilterLossesNoProximity(std::map<std::string, std::string> models);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterDimensions : public MagneticFilter {
    public:
        MagneticFilterDimensions() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

/**
 * @class MagneticFilterTurnCount
 * @brief Manufacturability / copper-burden proxy: turns × core size.
 *
 * Score = (Σ N_i) × core_width   (lower is better, caller inverts).
 *
 * The core-size factor is mandatory. For a fixed |Z| target on a CMC,
 * Z = µ_complex · ω · N² · Aₑ / lₑ, so N ∝ 1/√(µ·Aₑ/lₑ): bigger cores
 * minimise N. Scoring N alone therefore monotonically rewards the
 * largest core in the catalogue — which is the opposite of "manufacturable"
 * (more wire per turn, more copper, more cost). Multiplying by
 * core_width (OD for toroids, A-dim for two-piece sets) gives a simple
 * 1-D wire-length proxy that lets a small T 25 with N=20 outrank a big
 * T 140 with N=13. See TestTopologyCmc::Test_Cmc_AdviserMustNotPickOversizedToroid_WizardDefaults.
 *
 * Read-only: requires N to be already populated on the coil's
 * functional_description (e.g. via add_initial_turns_by_inductance /
 * filterMinimumImpedance). Returns valid=true with score=0 if no turns
 * have been assigned yet (caller-side gate ensures ordering).
 */
class MagneticFilterTurnCount : public MagneticFilter {
    public:
        MagneticFilterTurnCount() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterCoreMinimumImpedance : public MagneticFilter {
    private:
        // Candidate GATING runs the fast (OneLayer) capacitance path, explicitly, for exactly
        // the reason 2047e169 gave for MagneticFilterImpedance -- which is a DIFFERENT class,
        // and was the only one that commit fixed. MKF d424c32e made the full energy-based
        // capacitance model Impedance's default, which is right for analysing ONE magnetic but
        // wrong here: build_magnetizing_tank WINDS THE COIL and runs the per-turn energy sum on
        // every call, and this filter calls calculate_impedance once per requirement frequency
        // per candidate, plus up to five Newton re-bumps, for every core in the catalogue.
        //
        // Measured on Test_CoreAdviser_DMC_Default_Wizard_Hang_Repro (10 s budget), with the
        // default-constructed member below: filterMinimumImpedance alone took 1,038,396 ms over
        // 4,834 candidates -- 17.3 of the run's 17.4 minutes; every other stage in the pipeline
        // summed to about 7 s. That is the 28-minute regression 2047e169 describes, still live,
        // because the DMC/CMC suppression pipeline gates on THIS filter, not on that one.
        //
        // Rank the catalogue with the fast model; analyse the chosen design with the full one.
        Impedance _impedanceModel{/*fastCapacitance=*/true};
    public:
        MagneticFilterCoreMinimumImpedance() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterAreaNoParallels : public MagneticFilter {
    private:
        int _maximumNumberParallels = defaults.maximumNumberParallels;
    public:
        MagneticFilterAreaNoParallels() {};
        MagneticFilterAreaNoParallels(int maximumNumberParallels);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        std::pair<bool, double> evaluate_magnetic(Winding winding, Section section);
};

/**
 * @brief ABT #1177 (WP8, DFM rule R1): scores a wire candidate by the parity of the layer
 * count its winding would land on inside the section.
 *
 * An odd layer count forces a drag-back - a bump, window loss, extra leakage, a 45 to 90
 * degree wire crossing some safety standards forbid, and manual tape work. A single layer
 * cannot drag back and is never penalised. The scoring is the penalty (0 for an acceptable
 * candidate, 1 for an odd one), so it is inverted like every other cost-shaped filter. The
 * filter never rejects a candidate: an odd layer count is manufacturable, just worse.
 *
 * The numbers behind the rule live in src/data/dfm_rules.json; the only one this filter
 * needs is "a single layer is exempt", which is structural rather than tunable.
 */
class MagneticFilterLayerParity : public MagneticFilter {
    public:
        MagneticFilterLayerParity() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        std::pair<bool, double> evaluate_magnetic(Winding winding, Section section);
        /// The number of layers the winding needs inside the section, or nullopt when the
        /// section or the wire does not carry the dimensions to work it out.
        static std::optional<size_t> calculate_number_layers(Winding winding, Section section);
};

class MagneticFilterAreaWithParallels : public MagneticFilter {
    public:
        MagneticFilterAreaWithParallels() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        std::pair<bool, double> evaluate_magnetic(Winding winding, Section section, double numberSections, double sectionArea, bool allowNotFit);
};

class MagneticFilterEffectiveResistance : public MagneticFilter {
    private:
        double _maximumCurrent;
    public:
        MagneticFilterEffectiveResistance() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        std::pair<bool, double> evaluate_magnetic(Winding winding, double effectivefrequency, double temperature);
};

class MagneticFilterProximityFactor : public MagneticFilter {
    private:
        double _maximumCurrent;
    public:
        MagneticFilterProximityFactor() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        std::pair<bool, double> evaluate_magnetic(Winding winding, double effectiveSkinDepth, double temperature);
};

/**
 * @class MagneticFilterWindability
 *
 * Can each winding's wire actually be bent around the corner of the former it is wound on?
 * A turn wraps the column exactly the way the flexibility test wraps its mandrel, so the
 * standards give the answer directly (see WireBend): the former's corner radius must be at
 * least the mandrel's. Two thresholds, both citations rather than tuning knobs --
 * IEC 60317-0-1 Table 6 (or -0-2 Table 6) is the bend the insulation must survive AT ALL, and
 * Table 7 (or -0-2 clause 9) the bend it must survive and then be HEAT SHOCKED, which is what a
 * coil gets when it is soldered, varnish-baked and cycled.
 *
 * The bend judged is the one the coil actually lays -- the turn tangent to the former, so
 * former corner + standoff. A candidate is INVALID below the flexibility floor: that wire cannot
 * be wound on that former without cracking its enamel, and the only way it could be is by
 * standing off the former, which is not the layup the layout assumed. Below the heat-shock floor
 * it stays valid but scores worse, in proportion to how far short the bend falls.
 *
 * Wire types the standards do not cover (litz, foil, planar) are passed through untouched rather
 * than judged by a rule that was not written for them.
 */
class MagneticFilterWindability : public MagneticFilter {
    public:
        MagneticFilterWindability() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterSolidInsulationRequirements : public MagneticFilter {
    private:
        double _maximumCurrent;
    public:
        MagneticFilterSolidInsulationRequirements() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        std::pair<bool, double> evaluate_magnetic(Winding winding, WireSolidInsulationRequirements wireSolidInsulationRequirements);
};

class MagneticFilterTurnsRatios : public MagneticFilter {
    public:
        MagneticFilterTurnsRatios() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterMaximumDimensions : public MagneticFilter {
    public:
        MagneticFilterMaximumDimensions() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        bool applies_to(Magnetic* magnetic) const override;
};

class MagneticFilterSaturation : public MagneticFilter {
    public:
        MagneticFilterSaturation() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterDcCurrentDensity : public MagneticFilter {
    public:
        MagneticFilterDcCurrentDensity() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterEffectiveCurrentDensity : public MagneticFilter {
    public:
        MagneticFilterEffectiveCurrentDensity() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

/**
 * @class MagneticFilterDatasheetLimits
 * @brief Gate a catalogue part against its OWN datasheet-published electrical
 *        limits (rated current / voltage / saturation current).
 *
 * Applies only to a part whose datasheet publishes at least one of those
 * limits (applies_to): designed magnetics, and catalogue parts stating only an
 * inductance, are not judged by it and get no score from it.
 *
 * Checks each published limit against the operating point, skipping limits
 * that are absent:
 *   - rated current vs winding current RMS: `ratedCurrents`, or when only the
 *     ΔT-qualified `ratedCurrentPoints` is given, its smallest current;
 *   - saturation current vs the largest winding peak: the smallest of
 *     `saturationCurrentPeak` and every `saturationCurrents` criterion;
 *   - rated AC voltage vs voltage RMS, rated DC voltage vs |voltage offset|.
 *
 * No analytical physics: operating values are read from `Inputs`
 * (current/voltage processed RMS/peak/offset), datasheet values from the MAS
 * model. score is the UTILISATION -- the largest operating/limit ratio across
 * the checked limits, 1.0 at a limit. Lower is better, like every other
 * filter's score, so the usual invert=true ranks the part with the most
 * margin first. (Until 2026-09 this returned the smallest headroom, higher =
 * better, which invert=true turned upside down: the part closest to its
 * ratings ranked best.) valid=false when any limit is exceeded.
 *
 * `electrical` is a vector (one entry per connection configuration). The entry
 * whose `numberTurns` matches the candidate coil's turns is used; if none
 * matches (or `numberTurns` is unset, or the part has no coil) the most
 * conservative entry — smallest published rated current — is used rather
 * than guessing.
 *
 * Note (ABT #19): the MKF-pinned MAS uses `ratedCurrents` (array). Catalogues
 * still on the older `ratedCurrent` (scalar) schema deserialise to
 * get_rated_currents()==nullopt, so this filter does not apply to them until
 * the consumer's catalogue is migrated (asgard-side follow-up).
 */
class MagneticFilterDatasheetLimits : public MagneticFilter {
    public:
        MagneticFilterDatasheetLimits() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        bool applies_to(Magnetic* magnetic) const override;
};

class MagneticFilterImpedance : public MagneticFilter {
    public:
        MagneticFilterImpedance() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterMagnetizingInductance : public MagneticFilter {
    public:
        MagneticFilterMagnetizingInductance() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        bool applies_to(Magnetic* magnetic) const override;
};

/**
 * @class MagneticFilterLeakageInductance
 * @brief Filter for evaluating and scoring magnetic designs based on leakage inductance.
 *
 * For Common Mode Chokes (CMCs), lower leakage inductance indicates tighter coupling
 * between windings, which is essential for effective common-mode rejection.
 * The coupling coefficient k = 1 - (Lk/Lm) should be close to 1.
 *
 * For transformers requiring low leakage (e.g., LLC resonant converters),
 * this filter validates against the leakage_inductance design requirement.
 *
 * @note Returns leakage inductance in Henries. Lower values score better when inverted.
 */
enum class LeakageInductanceFilterMode {
    // Score Lk/Lm, lower is better (common-mode chokes). The historical behaviour, and the default.
    MINIMIZE_LEAKAGE_RATIO,
    // Score the distance of Lk (source winding 0 to winding i + 1) to
    // designRequirements.leakageInductance[i]; candidates outside the band are invalid (ABT #1176).
    TARGET
};

class MagneticFilterLeakageInductance : public MagneticFilter {
    private:
        LeakageInductanceFilterMode _mode = LeakageInductanceFilterMode::MINIMIZE_LEAKAGE_RATIO;
        std::pair<bool, double> evaluate_minimize_leakage_ratio(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs);
        std::pair<bool, double> evaluate_target(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs);
    public:
        MagneticFilterLeakageInductance() {};
        explicit MagneticFilterLeakageInductance(LeakageInductanceFilterMode mode) : _mode(mode) {};
        LeakageInductanceFilterMode get_mode() const { return _mode; }
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterSkinLossesDensity : public MagneticFilter {
    public:
        MagneticFilterSkinLossesDensity() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        std::pair<bool, double> evaluate_magnetic(Winding winding, SignalDescriptor current, double temperature);
};

class MagneticFilterFringingFactor : public MagneticFilter {
    private:
        std::map<std::string, std::string> _models;
        MagneticEnergy _magneticEnergy;
        double _requiredMagneticEnergy;
        double _fringingFactorLitmit = 1.2;
        std::shared_ptr<ReluctanceModel> _reluctanceModel;
        // Memoizes get_gapping_by_fringing_factor (a ≤100-iteration bisection)
        // keyed by core GEOMETRY + fringing limit. That solver overwrites the gap
        // length with its own search variable and reads only gap/column geometry
        // (fringing factor is purely geometric — vacuum permeability, gap area,
        // gap length), so its result is independent of material, turns, and the
        // incoming gap length. The ferrite standard-cores path expands one shape
        // into many materials with identical geometry, so without this the bisection
        // is recomputed identically per material. Instance-scoped (not static) so
        // there is no cross-run / cross-test shared state.
        std::map<std::string, double> _fringingFactorCache;

    public:
        MagneticFilterFringingFactor() {};
        MagneticFilterFringingFactor(Inputs inputs);
        MagneticFilterFringingFactor(Inputs inputs, std::map<std::string, std::string> models);
        void set_fringing_factor_limit(double limit) { _fringingFactorLitmit = limit; }
        double get_fringing_factor_limit() const { return _fringingFactorLitmit; }
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterVolume : public MagneticFilter {
    public:
        MagneticFilterVolume() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        bool applies_to(Magnetic* magnetic) const override;
};

class MagneticFilterArea : public MagneticFilter {
    public:
        MagneticFilterArea() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        bool applies_to(Magnetic* magnetic) const override;
};

class MagneticFilterHeight : public MagneticFilter {
    public:
        MagneticFilterHeight() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        bool applies_to(Magnetic* magnetic) const override;
};

class MagneticFilterTemperatureRise : public MagneticFilter {
    private:
        MagneticFilterLossesNoProximity _magneticFilterLossesNoProximity;
    public:
        MagneticFilterTemperatureRise() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterLossesTimesVolume : public MagneticFilter {
    public:
        MagneticFilterLossesTimesVolume() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterVolumeTimesTemperatureRise : public MagneticFilter {
    private:
        MagneticFilterTemperatureRise _magneticFilterTemperatureRise;
    public:
        MagneticFilterVolumeTimesTemperatureRise() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterLossesTimesVolumeTimesTemperatureRise : public MagneticFilter {
    private:
        MagneticFilterTemperatureRise _magneticFilterTemperatureRise;
    public:
        MagneticFilterLossesTimesVolumeTimesTemperatureRise() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterLossesNoProximityTimesVolume : public MagneticFilter {
    private:
        MagneticFilterLossesNoProximity _magneticFilterLossesNoProximity;
    public:
        MagneticFilterLossesNoProximityTimesVolume() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterLossesNoProximityTimesVolumeTimesTemperatureRise : public MagneticFilter {
    private:
        MagneticFilterTemperatureRise _magneticFilterTemperatureRise;
        MagneticFilterLossesNoProximity _magneticFilterLossesNoProximity;
    public:
        MagneticFilterLossesNoProximityTimesVolumeTimesTemperatureRise() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagnetomotiveForce : public MagneticFilter {
    public:
        MagnetomotiveForce() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};


// // Nice to have in the future
// class MagneticFilterMaximumWeight : public MagneticFilter {
//     public:
//         MagneticFilterMaximumWeight() {};
//         std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
// };

class MagneticFilterTemperature : public MagneticFilter {
    double _maximumTemperature = 130.0;
    // Orchestrator, not a fixed model: selects per material from its available
    // volumetric-losses methods (Steinmetz family, proprietary, loss factor)
    CoreLosses _coreLosses;
    MagnetizingInductance _magnetizingInductance;
public:
    MagneticFilterTemperature() {};
    MagneticFilterTemperature(Inputs inputs, double maximumTemperature);
    std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs,
                                              std::vector<Outputs>* outputs = nullptr);
};

// // Nice to have in the future
// class MagneticFilterStrayCapacitance : public MagneticFilter {
//     public:
//         MagneticFilterStrayCapacitance() {};
//         std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
// };

} // namespace OpenMagnetics