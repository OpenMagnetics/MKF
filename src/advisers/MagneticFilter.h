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

// Helpers the datasheet filters share for judging a datasheet coupled inductor by its ampere-turn
// current (MagneticFilterDatasheet.cpp).
std::vector<double> design_turns_ratios(const Inputs& inputs);
std::vector<IsolationSide> design_isolation_sides(const Inputs& inputs, Magnetic* magnetic);

// Key of the per-advise inductance/flux cache (cache_inductance_flux): the reference plus what
// the inductance and flux depend on that a candidate can change without changing its reference
// (turns of every winding, every gap length, the core material). Keyed on the reference alone, a
// candidate whose turns or gap an adviser step moved (the loss-optimal turns) read the flux of
// its earlier design.
std::string inductance_flux_cache_key(const Magnetic& magnetic);

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

        // Why the core losses of `magnetic` cannot be evaluated at all, or nullopt when they can: its
        // core material carries no core-loss data any model can run (e.g. `volumetricLosses: {}`).
        // Such a part is outside what every filter built on core losses can say -- they do not apply
        // to it (MagneticFilterCoreLossesBased) -- and the catalogue adviser flags it with this
        // reason instead of aborting the search on it or dropping it. Nullopt for a part without a
        // core: that is the datasheet-only case, which applies_to already answers.
        static std::optional<std::string> core_losses_not_evaluable_reason(Magnetic* magnetic);

        // Whether evaluating this filter computes core losses at the operating points. The catalogue
        // adviser applies its loss-model frequency-span gate (MagneticFilterLossModelFrequencySpan)
        // only when some filter of the flow, or its final simulation, computes them: a flow that
        // computes no losses has no reason to exclude a part for its loss fit (ABT #1679).
        virtual bool computes_core_losses() const { return false; }

        // The requirement frequencies this filter judged `magnetic` on, or nullopt for a filter
        // whose requirement is not frequency-wise. A filter that skips the points its data cannot
        // evaluate (MagneticFilterImpedance outside the material's mu(f) range) lists only the
        // points it judged, so the caller can say what the ranking rests on.
        virtual std::optional<std::vector<double>> get_judged_frequencies(Magnetic* magnetic, Inputs* inputs) const { return std::nullopt; }

        // The subset of get_judged_frequencies judged on the part's own measured data instead of the
        // model (MagneticFilterImpedance: points outside the material's mu(f) range judged on the
        // datasheet common-mode |Z|), or nullopt for a filter whose requirement is not frequency-wise.
        virtual std::optional<std::vector<double>> get_measured_frequencies(Magnetic* magnetic, Inputs* inputs) const { return std::nullopt; }

        // A copy of `magnetic` whose coil carries every winding the inputs excite. At the core
        // stage a multi-winding candidate carries a one-winding stand-in coil; the windings it
        // lacks are added exactly as the core adviser completes its results (correct_windings:
        // turns from the turns ratios, each winding's parallels sized to its own current). A
        // coil that already has every winding is returned as is, except a stand-in (its bobbin
        // still a name): the common-mode choke stand-in carries every winding, but only the first
        // has its seeded turns (the others hold the placeholder 1), so it is completed the same
        // way. One with more windings than the inputs excite throws.
        static Magnetic with_every_winding(const Magnetic& magnetic, const Inputs& inputs);

        // Conducting area of the round conductor whose own effective current density for
        // `current` is maximumEffectiveCurrentDensity: log bisection between rms / maximum (the
        // effective density is never below the DC one) and upperArea (a conductor known to be
        // enough). Throws when the current has no processed rms.
        static double get_conducting_area_for_current(const SignalDescriptor& current, double upperArea, double temperature,
                                                      double maximumEffectiveCurrentDensity);

        // A copy of a stand-in `magnetic` whose single-strand windings carry their current with
        // copper to spare (one stand-in strand, two skin depths: 18.6 mm at 50 Hz, below the
        // maximum effective current density at every operating point) are wound with a round
        // conductor the coil stage could pick: at least the copper MagneticFilterWindowCopperCapacity
        // counts (get_conducting_area_for_current), at most the strand, and as thick as the window
        // holds (filled at the capacity screen's utilisation, then the largest that fast_wind lays
        // out as whole layers). The core-stage losses and temperature then lay out and score copper
        // the coil stage can wind, not a strand no window can hold, nor the thinnest wire allowed
        // (a 32 A PFC choke at 12 A/mm2 dissipates ~200 W in its copper). Windings needing more
        // than one strand get more parallels of it by the same window-filling scale. With
        // fillWindow false only the spare strands are replaced, by the minimum copper (no layout
        // is tried). The coil must carry every winding the inputs excite.
        static Magnetic with_copper_sized_to_current(const Magnetic& magnetic, const Inputs& inputs, double maximumEffectiveCurrentDensity,
                                                     bool fillWindow = true);

        // A copy of `magnetic` whose windings' parallel strands are merged into one round
        // conductor of their copper area (same turns, same copper, one bundle): the coil the
        // core-stage thermal network is solved on, one node per merged turn.
        static Magnetic with_merged_strands(const Magnetic& magnetic, double temperature);
};

// Completes each candidate's one-winding stand-in coil with the windings the inputs excite
// (turns from the requirement turns ratios, parallels sized to each winding's current).
void correct_windings(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, const Inputs& inputs);

// A filter whose verdict is built on the core losses. It applies to what MagneticFilter does, and
// only when the core material has a core-loss model: for a part whose material has none it records
// no score (never a stand-in loss), and the part is ranked on the filters that do apply.
class MagneticFilterCoreLossesBased : public MagneticFilter {
    public:
        bool applies_to(Magnetic* magnetic) const override {
            return MagneticFilter::applies_to(magnetic) && !core_losses_not_evaluable_reason(magnetic);
        }
        bool computes_core_losses() const override { return true; }
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
        // The required area product reads the flux density from the material's loss model at each
        // operating frequency.
        bool computes_core_losses() const override { return true; }
        // The share of a core's winding window the bobbin leaves to the winding: the bobbin
        // model's factor for bobbin-wound cores, the inner-to-outer circumference rule for
        // toroids, 1 for printed windings.
        static double get_bobbin_filling_factor(const Core& core, std::optional<WiringTechnology> wiringTechnology);
};

/**
 * @class MagneticFilterWindowCopperCapacity
 * @brief Core-stage copper screen: can the winding window hold the copper of every winding?
 *
 * At the core stage a candidate carries a one-winding stand-in coil; the other windings
 * follow from the turns ratios. Each winding needs, per turn, the strands of the stand-in's
 * wire (two skin depths, the finest round strand the wire adviser's litz is built from)
 * that keep its effective current density at or below the maximum
 * (the effective current density of one strand over the maximum, rounded up). When one
 * strand is more than enough (always at line frequency, where two skin depths are
 * centimetres) the coil stage can pick a thinner wire, so only that fraction of a strand is
 * counted. The window holds that copper at the utilisation MagneticFilterAreaProduct sizes
 * cores with: the round-wire filling factor of each winding's conductor times the bobbin
 * filling factor. A core whose
 * window cannot hold all windings' copper cannot be wound by the coil stage within its
 * current-density limit, whatever wire it picks, and is rejected here with its numbers.
 * This is a necessary condition only (litz of the finest strands at the bobbin's best fill):
 * it never rejects a core the coil stage could wind.
 * Printed windings are not screened: their copper is the PCB's, not strands in a window.
 */
class MagneticFilterWindowCopperCapacity : public MagneticFilter {
    private:
        std::vector<double> _turnsRatios;
        double _temperature = 0;
        double _maximumEffectiveCurrentDensity = 0;
        std::string _lastReason;
        // Copper area of the conductor sized to one winding's current, by (strand diameter, winding, operating point).
        std::map<std::tuple<double, size_t, size_t>, double> _conductorAreaCache;
    public:
        MagneticFilterWindowCopperCapacity() {};
        MagneticFilterWindowCopperCapacity(Inputs inputs, double maximumEffectiveCurrentDensity);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        // Why the last evaluated candidate was rejected (empty when it passed).
        const std::string& get_last_reason() const { return _lastReason; }
        static bool applies_to(const Inputs& inputs);
};

/**
 * @class MagneticFilterWireWithinLimits
 * @brief Can every winding of a completed coil get a wire within the coil stage's limits?
 *
 * The window copper capacity is necessary, not sufficient: the coil stage gives each winding a
 * section of the window (CoilAdviser::get_advised_sections) and looks for a wire that carries
 * the winding's current at or below the maximum effective current density with at most the
 * maximum number of parallels and fits that section (WireAdviser, with litz synthesis, as the
 * coil stage runs it). A 2-turn 33 A secondary in half an E 20 window found no such wire, so
 * every coil the coil stage made for that core came back INVALID. This filter runs that same
 * search for each winding on each pattern the coil stage would try, at its fewest repetitions
 * and without margin tape (the roomiest sections): a core where no pattern gives every winding
 * a wire is rejected, naming the winding. The candidate's coil must carry every winding.
 */
class CoilAdviser;
class WireAdviser;
class MagneticFilterWireWithinLimits : public MagneticFilter {
    private:
        double _maximumEffectiveCurrentDensity = 0;
        int _maximumNumberParallels = 0;
        std::string _lastReason;
        std::shared_ptr<CoilAdviser> _coilAdviser;
        std::shared_ptr<WireAdviser> _wireAdviser;
        // The wires the coil stage advises from (CoilAdviser::get_catalogue_wires), built once.
        std::optional<std::vector<Wire>> _catalogueWires;
    public:
        MagneticFilterWireWithinLimits(double maximumEffectiveCurrentDensity, int maximumNumberParallels);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        const std::string& get_last_reason() const { return _lastReason; }
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

// Unit-cost basis of the COST filter: US$ per piece at a 1000-piece production quantity. Every value
// carries the source it was read from; a quantity with no source is marked as an assumption.
struct CostBasis {
    // Refined copper, LME cash settlement US$14,737/t on 2026-09-08 (Cochilco, reported by Reporte Minero,
    // https://www.reporteminero.cl/noticia/noticias/2026/09/precio-cobre-nuevo-record-alza-lme-septiembre-2026).
    // The magnet-wire conversion premium over refined copper has no source yet and is not added.
    double copperPricePerKg = 14.737;
    // Litz strand premium per unit copper mass, C_m(d) = 1 + k1/d^6 + k2/d^2: Sullivan, "Cost-constrained
    // selection of strand diameter and number in a litz-wire transformer winding", IEEE TPEL 16(2) 2001,
    // eq. 4, with the constants refit to newer prices in Sullivan & Zhang, "Simplified design method for
    // litz wire", APEC 2014 (https://www.ryz.web.illinois.edu/pdf/simplitz_apec2014.pdf).
    double litzK1 = 6e-26;   // m^6
    double litzK2 = 2.7e-9;  // m^2
    // Labour rate: Mexico, semi-skilled, fully loaded, US$7.27/h (Tetakawi, 2026, a shelter-services vendor).
    double laborRatePerHour = 7.27;
    // Winding time: 1.05 min for a 2 x 5-turn machine-assisted toroid (Kamil et al., IJIM 2021), i.e.
    // 0.105 min per turn. ASSUMPTION (unsourced): the time scales linearly with the number of turns and
    // applies to bobbin winding too.
    double windingMinutesPerTurn = 0.105;
    // Winding is 50 % of the labour (termination 40 %, encapsulation 10 %): US patent 5,781,091 (1995).
    double windingShareOfLabour = 0.5;
    // Core price per set (or per toroid) as a power law of the core mass, price = c * mass^k (US$, kg),
    // one law per MAS material type and materialComposition, keyed "type/composition" (e.g. "ferrite/MnZn",
    // "powder/FeSiAl"; a "proprietary" composition adds the manufacturer, e.g. "powder/proprietary/Micrometals",
    // since the word names no alloy). Price does not
    // scale with mass: small cores cost several times more per kg, so a flat US$/kg misprices them. The
    // values are fitted where they are set. A composition with no fitted law throws.
    struct PowerLaw {
        double coefficient;  // c, US$ for a 1 kg core
        double exponent;     // k
    };
    std::map<std::string, PowerLaw> corePrice;
};

struct CostBreakdown {
    // Empty when the core has no price law (its material has no composition, or no law is fitted for its
    // type and composition); coreUnpricedReason then says why. Never estimated.
    std::optional<double> core;
    double conductor;
    double labour;
    std::string coreUnpricedReason;
    bool priced() const { return core.has_value(); }
    // The unit cost; throws for an unpriced core, whose cost is unknown.
    double total() const;
};

class MagneticFilterCost : public MagneticFilter {
    public:
        MagneticFilterCost() {};
        explicit MagneticFilterCost(CostBasis basis) : _basis(std::move(basis)) {};
        // Scores the unit cost in US$. A core with no price law is not dropped and not estimated: it returns
        // {false, 0}, which a non-strict adviser flow ranks worst for cost (the ABT #801 failed-scoring rule).
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        CostBreakdown calculate_cost(Magnetic& magnetic) const;
        static CostBasis default_basis();
    private:
        std::optional<CostBasis> _basis;
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

class MagneticFilterCoreAndDcLosses : public MagneticFilterCoreLossesBased {
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

class MagneticFilterCoreDcAndSkinLosses : public MagneticFilterCoreLossesBased {
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
        // A multi-winding candidate still carrying its one-winding stand-in coil: core losses plus
        // the DC and skin-effect losses of EVERY winding, on the coil completed from the turns
        // ratios (with_every_winding) and fast-wound. No turns sweep (the ratios fix the turns)
        // and the candidate keeps its stand-in coil.
        std::pair<bool, double> evaluate_stand_in_with_every_winding(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterLosses : public MagneticFilterCoreLossesBased {
    private:
        std::map<std::string, std::string> _models;
        MagneticSimulator _magneticSimulator;
    public:
        MagneticFilterLosses() {};
        MagneticFilterLosses(std::map<std::string, std::string> models);
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterLossesNoProximity : public MagneticFilterCoreLossesBased {
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
        // The same verdict for a bare wire, without building a Winding around it. Scores nothing:
        // every candidate this filter keeps scores 0, so it is a pure per-wire predicate.
        bool wire_fits(Wire& wire, int64_t numberParallels, int64_t numberTurns, const Section& section) const;
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
    private:
        // ABT #1699: wires the real winder refused in the strict pass (real winding on only).
        size_t _realWindingRefusals = 0;
        std::string _firstRealWindingRefusal;
    public:
        MagneticFilterAreaWithParallels() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        std::pair<bool, double> evaluate_magnetic(Winding winding, Section section, double numberSections, double sectionArea, bool allowNotFit);
        size_t get_real_winding_refusals() const { return _realWindingRefusals; }
        const std::string& get_first_real_winding_refusal() const { return _firstRealWindingRefusal; }
};

/**
 * @class MagneticFilterRealWinding
 * @brief Does the real winder build this wound magnetic as it is stored? (ABT #1699)
 *
 * A coil an adviser returns is stored and rebuilt downstream from its JSON with real winding
 * geometry (the web's 3D view, MVB++ magnetic_autocomplete_safe). This filter asks the real
 * winder itself, through real_winding_refusal: the magnetic rebuilt from its stored form (its
 * functional and section description; layers and turns are laid out again) must
 * wind, apply its connection blocking and keep every turn inside the winding window. The
 * refusal the real winder gives is kept as the last reason. Planar coils have no real-winding
 * model (ABT #492): evaluating one throws.
 */
class MagneticFilterRealWinding : public MagneticFilter {
    private:
        std::string _lastReason;
    public:
        MagneticFilterRealWinding() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        // Why the last evaluated magnetic was refused (empty when it passed).
        const std::string& get_last_reason() const { return _lastReason; }
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

/**
 * @class MagneticFilterLossModelFrequencySpan
 *
 * Gate, not a score: can the core loss model MKF would use for this material be evaluated at every
 * operating frequency? The Steinmetz ranges in MAS are fits over [minimumFrequency, maximumFrequency];
 * outside that span CoreLossesModel::get_steinmetz_coefficients throws MaterialFrequencyOutOfSpanException
 * rather than extrapolate (ABT #1456: nine 0.5-5 MHz grades read at 100 kHz gave 2-18 % of 3C95's loss and
 * won every ranking). Advisers must drop such materials through this filter, which returns {false, 0},
 * instead of letting the loss filters throw: the MagneticAdviser aborts after 8 identical throws.
 *
 * Only the Steinmetz family (Steinmetz, iGSE, ciGSE, Barg, Albach, MSE, NSE) has a declared span. The other
 * methods (Magnetics / Micrometals / Poco / TDG closed forms, Roshen, loss factor) declare none in the
 * schema, so a material evaluated with one of them always passes.
 *
 * With Settings::allowMaterialDataExtrapolation on (explicit opt-in, off by default) the span does not
 * gate: the losses are evaluated outside it and every such evaluation logs a WARNING. Never under an
 * adviser (ABT #1652): while any adviser object exists the flag reads false and the span gates as usual.
 */
class MagneticFilterLossModelFrequencySpan : public MagneticFilter {
    public:
        MagneticFilterLossModelFrequencySpan() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        // Needs only the core's material, so it judges parts without a coil too. A material with no
        // core-loss model at all has no span to judge: the loss filters do not apply to it
        // (core_losses_not_evaluable_reason) and the catalogue adviser flags the part, rather than
        // this gate dropping it.
        bool applies_to(Magnetic* magnetic) const override { return magnetic->has_core() && !core_losses_not_evaluable_reason(magnetic); }
        // The same verdict for a bare material. `model` is the caller's requested core loss model; like
        // CoreLosses, it heads the settings' model order and the first model the material supports is
        // the one judged. A material no model can evaluate is not evaluable (false). When the Settings
        // carry an explicitly requested model (ABT #1497), a material that cannot run that model, or
        // not at every operating frequency, is not evaluable either.
        static bool is_material_evaluable(const CoreMaterial& material, const Inputs& inputs,
                                          std::optional<CoreLossesModels> model = std::nullopt);
        // Why the material is not evaluable (the model, its fitted span and the operating frequencies
        // outside it), or nullopt when it is. is_material_evaluable is its negation.
        static std::optional<std::string> material_not_evaluable_reason(const CoreMaterial& material, const Inputs& inputs,
                                                                        std::optional<CoreLossesModels> model = std::nullopt);
};

class MagneticFilterSolidInsulationRequirements : public MagneticFilter {
    private:
        double _maximumCurrent;
    public:
        MagneticFilterSolidInsulationRequirements() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        std::pair<bool, double> evaluate_magnetic(Winding winding, WireSolidInsulationRequirements wireSolidInsulationRequirements);
        // The same verdict and scoring for a bare wire, without building a Winding around it.
        std::pair<bool, double> evaluate_wire(Wire& wire, const WireSolidInsulationRequirements& wireSolidInsulationRequirements);
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
        // Design-mode checks: each axis as stated, no rotation, so a height limit stays a
        // height limit (the catalogue check above allows rotation). core_fits is a necessary
        // condition before any coil exists; magnetic_fits checks the wound assembly.
        static bool core_fits(Core& core, const Inputs& inputs);
        static bool magnetic_fits(Magnetic& magnetic, const Inputs& inputs);
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
        // Measured common-mode |Z| of a part, (frequency Hz, |Z| Ohm), strictly increasing in frequency.
        using MeasuredImpedanceCurve = std::vector<std::pair<double, double>>;

        MagneticFilterImpedance() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
        // The minimumImpedance frequencies the filter judges: those inside the core material's
        // tabulated mu(f) range (judged on the model), and those outside it but inside the part's
        // measured common-mode |Z| curve (judged on the measurement). Nullopt without a
        // minimumImpedance requirement.
        std::optional<std::vector<double>> get_judged_frequencies(Magnetic* magnetic, Inputs* inputs) const override;
        // The judged frequencies judged on the measured curve (outside the material's mu(f) range).
        std::optional<std::vector<double>> get_measured_frequencies(Magnetic* magnetic, Inputs* inputs) const override;

        // The part's zero-bias common-mode |Z| from its datasheet: the impedancePoints of the
        // manufacturerInfo.datasheetInfo.electrical entry with subtype commonModeChoke, keeping only
        // the points without a DC-bias current (or with current 0). Nullopt when the part carries
        // none. Throws InvalidDatasheetImpedanceException on a malformed point, two values at one
        // frequency, or more than one zero-bias curve (several entries, or per-winding curves).
        static std::optional<MeasuredImpedanceCurve> get_measured_impedance_curve(const Magnetic& magnetic);
        // |Z| at `frequency` from the measured curve: linear interpolation of log|Z| against log f
        // between the two bracketing measured points (the point itself when measured there). Throws
        // outside the measured frequency range: the curve is never extrapolated.
        static double interpolate_measured_impedance(const MeasuredImpedanceCurve& curve, double frequency);

    private:
        // Per requirement point: judged on the model (inside the material's mu(f) range), on the
        // measured curve (outside it, inside the curve's range), or not judged.
        enum class ImpedanceSource { MODEL, MEASURED, NOT_JUDGED };
        struct RequirementCoverage {
            std::vector<ImpedanceSource> sources;  // one per minimumImpedance point, in order
            std::optional<MeasuredImpedanceCurve> measuredCurve;
            double minimumMaterialFrequency;
            double maximumMaterialFrequency;
        };
        static RequirementCoverage classify_requirement(Magnetic* magnetic, Inputs* inputs);
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

/**
 * @class MagneticFilterInductorTurnsAndGapByLosses
 * @brief ABT #1426: moves a gapped inductor candidate from its saturation-floor (N, gap) to
 * the pair that minimises the estimated core + copper loss, never below the floor.
 *
 * Not a gate: evaluate_magnetic always returns valid and rewrites the candidate's turns and
 * gap in place when a buildable lower-loss pair exists (see MagneticFilterLosses.cpp for
 * the model and the search). The returned score is the estimated mean total loss (W) of
 * the pair it settled on, 0 when the candidate is outside its scope.
 */
class MagneticFilterInductorTurnsAndGapByLosses : public MagneticFilterCoreLossesBased {
    private:
        std::map<std::string, std::string> _models;
        MagnetizingInductance _magnetizingInductance;
        std::shared_ptr<CoreLossesModel> _coreLossesModelSteinmetz;
        std::shared_ptr<CoreLossesModel> _coreLossesModelProprietary;
        MagneticFilterFringingFactor _fringingFactorFilter;
        MagneticFilterMagnetizingInductance _inductanceFilter;
        MagneticFilterSaturation _saturationFilter;
        double _targetInductance = 0;
        std::vector<double> _temperatures;
        std::vector<double> _primaryReferredCurrentsRms;
        std::vector<std::vector<double>> _primaryReferredCurrentHarmonicAmplitudes;
        std::vector<std::vector<double>> _primaryReferredCurrentHarmonicFrequencies;
        double _maximumEffectiveFrequency = 0;
        std::vector<double> _turnsRatios;
        std::vector<std::vector<SignalDescriptor>> _windingCurrents;

        // Mean over the operating points of the core losses at the target inductance and
        // numberTurns (the flux L i / N the re-solved gap delivers). NaN when the core has
        // no loss model to evaluate (neither Steinmetz nor a proprietary formula).
        double calculate_core_losses(const Core& core, double numberTurns, std::vector<OperatingPoint>& preparedOperatingPoints);

    public:
        MagneticFilterInductorTurnsAndGapByLosses(Inputs inputs, std::map<std::string, std::string> models);
        // Designs this applies to: energy-storing inductors with an inductance target, not
        // suppression chokes (sized by impedance).
        static bool applies_to_design(const Inputs& inputs);
        // Candidates this applies to: a discrete-gap core (not a toroid, not PQI/UI) with the
        // adviser's one-winding stand-in coil.
        static bool applies_to_candidate(const Magnetic& magnetic);
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

class MagneticFilterTemperatureRise : public MagneticFilterCoreLossesBased {
    private:
        MagneticFilterLossesNoProximity _magneticFilterLossesNoProximity;
    public:
        MagneticFilterTemperatureRise() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterLossesTimesVolume : public MagneticFilterCoreLossesBased {
    public:
        MagneticFilterLossesTimesVolume() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterVolumeTimesTemperatureRise : public MagneticFilterCoreLossesBased {
    private:
        MagneticFilterTemperatureRise _magneticFilterTemperatureRise;
    public:
        MagneticFilterVolumeTimesTemperatureRise() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterLossesTimesVolumeTimesTemperatureRise : public MagneticFilterCoreLossesBased {
    private:
        MagneticFilterTemperatureRise _magneticFilterTemperatureRise;
    public:
        MagneticFilterLossesTimesVolumeTimesTemperatureRise() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterLossesNoProximityTimesVolume : public MagneticFilterCoreLossesBased {
    private:
        MagneticFilterLossesNoProximity _magneticFilterLossesNoProximity;
    public:
        MagneticFilterLossesNoProximityTimesVolume() {};
        std::pair<bool, double> evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs = nullptr);
};

class MagneticFilterLossesNoProximityTimesVolumeTimesTemperatureRise : public MagneticFilterCoreLossesBased {
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

/**
 * @class MagneticFilterTemperature
 * @brief Core-stage temperature gate: the ThermalNetwork's hottest point at each operating point.
 *
 * The solve carries the core losses AND the copper: the DC and skin-effect losses of every
 * winding, on the candidate's coil completed from the turns ratios (with_every_winding) and
 * fast-wound, per turn. It used to solve the core alone, so a transformer whose copper dominates
 * (a 10 kW PSFB) looked cool here and every candidate then failed hot after the 8-12 s coil
 * stage. Proximity losses are not included (they need the coil stage's real layout), so the
 * estimate is still a lower bound on the final temperature. PQI / UI shapes, whose integrated
 * windings fast_wind() does not lay out, are solved core-only (the loss filters' policy).
 */
class MagneticFilterTemperature : public MagneticFilterCoreLossesBased {
    double _maximumTemperature = 130.0;
    // Orchestrator, not a fixed model: selects per material from its available
    // volumetric-losses methods (Steinmetz family, proprietary, loss factor)
    CoreLosses _coreLosses;
    MagnetizingInductance _magnetizingInductance;
    WindingOhmicLosses _windingOhmicLosses;
    WindingSkinEffectLosses _windingSkinEffectLosses;
    bool _sizeStandInCopper = false;
public:
    MagneticFilterTemperature() {};
    MagneticFilterTemperature(Inputs inputs, double maximumTemperature);
    // The core adviser's candidates carry a placeholder stand-in coil: their copper is laid out
    // and scored as the coil stage could wind it (with_copper_sized_to_current). Off by default:
    // a real coil is solved with its own wires.
    void set_size_stand_in_copper(bool value) { _sizeStandInCopper = value; }
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