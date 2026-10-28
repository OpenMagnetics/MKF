#include "support/NumberFormat.h"
#include "physical_models/AmplitudePermeability.h"
#include "physical_models/InitialPermeability.h"
#include "physical_models/CurveBuildScope.h"
#include "support/Exceptions.h"
#include <magic_enum.hpp>
#include <algorithm>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <tuple>


namespace OpenMagnetics {

// ============================================================================
// Remanence and coercive force: MAS value, or a flagged class estimate (ABT #1951).
// Alf, 2026-10-18: "an approved class-based estimator, like the heat-conductivity estimator of
// #1752/#1753: Br from Bs and the material class's loop type (round/flat/square), Hc likewise from
// the class; sourced data always wins; throw when the class is unknown."
//
// Derivation (MAS da80482, 1073 core materials; values at the tabulated point nearest 25 C, Br/Bs
// against the material's own Bs at that temperature; one point per document that gives all its
// grades one value, else one per grade; medians):
//   pool                        Br/Bs   (points, documents)   Hc A/m (points, documents)   loop
//   class FERRITE MN_ZN         0.255   (297, 110)            12.7   (295, 109)            round
//   class FERRITE NI_ZN         0.600   (138, 96)             58     (181, 96)             round
//   class FERRITE MG_ZN         0.621   (6, 6)                25     (6, 6)                round
//   material type FERRITE       0.326   (441, 211)            16     (482, 210)            round
//   NANOCRYSTALLINE: split by loop family, see below.
//   class AMORPHOUS FE_SI       0.815   (1: Metglas 2605SA1)  -                            square
//   material type AMORPHOUS     0.728   (2, 2: AF, Metglas 2605SA1)  2.4 (3, 3)            square
//   material type POWDER        0.0082  (91, 84)              151    (90, 84)              flat (distributed gap)
//   powder classes, Br/Bs / Hc: FE_SI_AL 0.0068 / 79.6, FE_SI 0.0066 / 267, FE_NI 0.0145 / 191,
//   FE_NI_MO 0.0122 / 63.7, IRON 0.0202 / 398, CARBONYL_IRON 0.0028 / 239.
//   ELECTRICAL_STEEL: Hc 0.755 (7 grades, 1 document) but no Br anywhere: a steel without Br throws.
// The table is what the code below computes from MAS; it is not used as data.
//
// NANOCRYSTALLINE is split by loop family (Alf, 2026-10-19). A nanocrystalline grade's permeability is
// set by its field anneal, which also sets its loop: the low-permeability grades are transverse-field
// annealed, with a flat (sheared) loop and a remanence of a few mT (Nanoperm 80000: Br 0.31 T at Hc 3
// A/m; a Nanoperm 1000 loop through the same Hc holds about mu0 x 1000 x 3 A/m = 4 mT), while the
// high-permeability grades have a round loop near Br/Bs 0.5. A pooled Br/Bs of the round grades is
// therefore not an estimate for a flat-loop grade. MAS core materials carry no loop-shape field, so the
// family is read from the initial permeability the record tabulates (InitialPermeability at the default
// ambient temperature, no bias, no frequency: the mu_i column below). The 18 nanocrystalline grades
// of MAS c1b29c0, sorted by mu_i (Br and Hc: sourced values at 25 C):
//   grade              mu_i      Br sourced   Hc sourced   ratio to the next grade's mu_i
//   Nanoperm 1000        998     -            3 A/m        2.0
//   Nanoperm 2000       2004     -            3 A/m        2.0
//   Nanoperm 4000       4044     -            3 A/m        2.1
//   Nanoperm 8000       8293     -            3 A/m        3.5   <- the only ratio above 2.1
//   SC-1K107-LP        29342     -            -            1.0
//   Nanoperm 30000     29787     -            3 A/m        1.3
//   VITROPERM 500F     40000     -            3.0 A/m      1.0
//   VITROPERM 550 HF   40000     -            3.0 A/m      1.1
//   FT-3M              42685     0.60 T       0.8 A/m      1.4
//   Nanoperm 60000     60000     -            -            1.3
//   AN, 1K107, Yunlu YNB 80000, Nanoperm 80000 80153 (Br: AN 0.60 T, Nanoperm 80000 0.31 T; Hc 0.8,
//     1.5, 0.8, 3 A/m); Nanoperm 85000 85000; Nanoperm 90000 87339 (Hc 3 A/m); SC-1K107 87998;
//     Nanoperm 100000 100000. All ratios above Nanoperm 8000 are 1.4 or less.
// The boundary sits in that gap: FLAT_LOOP_LOW_PERMEABILITY below 15000 (Nanoperm 1000 to 8000),
// ROUND_LOOP_HIGH_PERMEABILITY from 15000 up. Each family pools only its own sourced members and never
// borrows from the other: no flat-loop grade has a sourced remanence, so the Nanoperm 1000/2000/4000/
// 8000 remanence throws; their coercive force is their own (3 A/m). Round family pools, Br/Bs 0.48 (3
// points, 3 documents: AN, FT-3M, Nanoperm 80000), Hc 1.15 A/m (6 points, 6 documents).
// ============================================================================
namespace {

enum class NanocrystallineLoop { FLAT_LOOP_LOW_PERMEABILITY, ROUND_LOOP_HIGH_PERMEABILITY };

// The initial permeability that separates the two nanocrystalline loop families (derivation above).
constexpr double nanocrystallineFlatLoopPermeabilityBoundary = 15000;

std::optional<NanocrystallineLoop> nanocrystalline_loop(const CoreMaterial& material) {
    if (material.get_material() != MaterialType::NANOCRYSTALLINE) {
        return std::nullopt;
    }
    const double initialPermeability = InitialPermeability::get_initial_permeability(material, Defaults().ambientTemperature);
    if (!std::isfinite(initialPermeability) || initialPermeability <= 0) {
        throw MaterialException(ErrorCode::MATERIAL_INVALID_PROPERTY,
            "Core material '" + material.get_name() + "' has an initial permeability of " + std::to_string(initialPermeability) +
            " at the default ambient temperature: its nanocrystalline loop family cannot be told.");
    }
    return initialPermeability < nanocrystallineFlatLoopPermeabilityBoundary ? NanocrystallineLoop::FLAT_LOOP_LOW_PERMEABILITY
                                                                             : NanocrystallineLoop::ROUND_LOOP_HIGH_PERMEABILITY;
}

std::string loop_family_label(const std::optional<NanocrystallineLoop>& family) {
    if (!family) {
        return "";
    }
    std::ostringstream label;
    label << (family.value() == NanocrystallineLoop::FLAT_LOOP_LOW_PERMEABILITY ? ", flat loop (initial permeability below "
                                                                                 : ", round loop (initial permeability from ")
          << nanocrystallineFlatLoopPermeabilityBoundary << ")";
    return label.str();
}

enum class LoopProperty { REMANENCE, COERCIVE_FORCE };

const char* loop_property_field(LoopProperty property) {
    return property == LoopProperty::REMANENCE ? "remanence" : "coerciveForce";
}

const char* loop_property_name(LoopProperty property) {
    return property == LoopProperty::REMANENCE ? "remanence" : "coercive force";
}

bool has_loop_property(const CoreMaterial& material, LoopProperty property) {
    return property == LoopProperty::REMANENCE ? material.get_remanence().has_value() : material.get_coercive_force().has_value();
}

// The documents a value is sourced from. A non-retracted provenance entry naming the field decides:
// its documents, unless it only says the value was derived (not a source). A value no entry names is
// the maker's catalogue value when its record cites the maker's datasheet.
std::vector<std::string> loop_property_sources(const CoreMaterial& material, LoopProperty property) {
    std::vector<std::string> sources;
    if (!has_loop_property(material, property)) {
        return sources;
    }
    const std::string field = loop_property_field(property);
    bool named = false;
    if (material.get_provenance()) {
        const auto provenance = material.get_provenance().value();
        for (const auto& entry : provenance) {
            if (entry.get_retracted() && entry.get_retracted().value()) {
                continue;
            }
            if (!entry.get_fields()) {
                continue;
            }
            const auto fields = entry.get_fields().value();
            if (std::find(fields.begin(), fields.end(), field) == fields.end()) {
                continue;
            }
            named = true;
            if (entry.get_source() == ProvenanceSource::DERIVED) {
                continue;
            }
            if (entry.get_source_url()) {
                sources.push_back(entry.get_source_url().value());
            }
            else if (entry.get_source_name()) {
                sources.push_back(entry.get_source_name().value());
            }
        }
    }
    if (!named && material.get_manufacturer_info().get_datasheet_url()) {
        sources.push_back(material.get_manufacturer_info().get_datasheet_url().value());
    }
    return sources;
}

// The member's pooled value: Br/Bs or Hc at its tabulated point nearest the default ambient
// temperature (ties: the lower temperature). Disengaged when a Br has no Bs to divide by.
std::optional<double> loop_property_pool_value(const CoreMaterial& material, LoopProperty property) {
    const double reference = Defaults().ambientTemperature;
    std::vector<std::pair<double, double>> points;  // (temperature, value)
    if (property == LoopProperty::REMANENCE) {
        const auto remanenceData = material.get_remanence().value();
        for (const auto& datum : remanenceData) {
            points.push_back({datum.get_temperature(), datum.get_magnetic_flux_density()});
        }
    }
    else {
        const auto coerciveForceData = material.get_coercive_force().value();
        for (const auto& datum : coerciveForceData) {
            points.push_back({datum.get_temperature(), datum.get_magnetic_field()});
        }
    }
    if (points.empty()) {
        return std::nullopt;
    }
    auto nearest = *std::min_element(points.begin(), points.end(), [reference](const auto& a, const auto& b) {
        const double da = std::abs(a.first - reference);
        const double db = std::abs(b.first - reference);
        return da < db || (da == db && a.first < b.first);
    });
    double value = nearest.second;
    if (property == LoopProperty::REMANENCE) {
        if (material.get_saturation().empty()) {
            return std::nullopt;
        }
        value /= Core::get_magnetic_flux_density_saturation(material, nearest.first, false);
    }
    if (!std::isfinite(value) || value <= 0) {
        std::ostringstream message;
        message << "Core material '" << material.get_name() << "' has a malformed " << loop_property_name(property)
                << " in MAS (pooled value " << value << (property == LoopProperty::REMANENCE ? " Br/Bs" : " A/m")
                << " at " << nearest.first << " C): it must be a finite positive number. Fix the record in MAS.";
        throw MaterialException(ErrorCode::MATERIAL_INVALID_PROPERTY, message.str());
    }
    return value;
}

struct LoopMember {
    std::string name;
    MaterialType materialType;
    std::optional<MaterialComposition> composition;
    std::vector<std::string> sources;
    double value;
    std::optional<NanocrystallineLoop> loopFamily;
};

struct LoopPool {
    double value;
    std::vector<std::string> sourceMaterials;
    size_t numberDocuments;
    size_t numberDataPoints;
};

double median(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    const size_t n = values.size();
    return n % 2 == 1 ? values[n / 2] : (values[n / 2 - 1] + values[n / 2]) / 2.0;
}

// As the heat-conductivity estimator's pool_data_points: each member counts under its citing document
// that covers the most members of the pool (ties: the smallest identifier); a document that gives all
// its grades one value is ONE point, one that gives each grade its own value is one point per grade.
std::optional<LoopPool> pool_median(const std::vector<const LoopMember*>& members) {
    if (members.empty()) {
        return std::nullopt;
    }
    std::map<std::string, size_t> citationsInPool;
    for (const auto* member : members) {
        for (const auto& source : std::set<std::string>(member->sources.begin(), member->sources.end())) {
            ++citationsInPool[source];
        }
    }
    std::map<std::string, std::vector<const LoopMember*>> membersBySource;
    for (const auto* member : members) {
        const std::string* widest = nullptr;
        for (const auto& source : member->sources) {
            if (!widest || citationsInPool.at(source) > citationsInPool.at(*widest) ||
                (citationsInPool.at(source) == citationsInPool.at(*widest) && source < *widest)) {
                widest = &source;
            }
        }
        membersBySource[*widest].push_back(member);
    }
    LoopPool pool;
    std::vector<double> points;
    for (const auto& [source, sourceMembers] : membersBySource) {
        std::set<double> distinct;
        for (const auto* member : sourceMembers) {
            distinct.insert(member->value);
            pool.sourceMaterials.push_back(member->name);
        }
        if (distinct.size() == 1) {
            points.push_back(sourceMembers.front()->value);
            continue;
        }
        for (const auto* member : sourceMembers) {
            points.push_back(member->value);
        }
    }
    std::sort(pool.sourceMaterials.begin(), pool.sourceMaterials.end());
    pool.numberDocuments = membersBySource.size();
    pool.numberDataPoints = points.size();
    pool.value = median(points);
    return pool;
}

// The sourced members of every material, per property; rebuilt when the material database changes size.
// generation: set to a number that changes every time a member list is rebuilt (the key of estimate_loop_property's memo).
const std::vector<LoopMember>& sourced_loop_members(LoopProperty property, size_t& generation) {
    static std::mutex membersMutex;
    static std::map<LoopProperty, std::pair<size_t, std::vector<LoopMember>>> cache;
    static size_t rebuilds = 0;
    std::lock_guard<std::mutex> lock(membersMutex);
    if (coreMaterialDatabase.empty()) {
        load_core_materials();
    }
    auto& entry = cache[property];
    if (entry.first == coreMaterialDatabase.size() && !entry.second.empty()) {
        generation = rebuilds;
        return entry.second;
    }
    ++rebuilds;
    std::vector<LoopMember> members;
    for (const auto& [materialName, material] : coreMaterialDatabase) {
        auto sources = loop_property_sources(material, property);
        if (sources.empty()) {
            continue;
        }
        const auto value = loop_property_pool_value(material, property);
        if (!value) {
            continue;
        }
        members.push_back({material.get_name(), material.get_material(), material.get_material_composition(), std::move(sources), value.value(),
                           nanocrystalline_loop(material)});
    }
    entry = {coreMaterialDatabase.size(), std::move(members)};
    generation = rebuilds;
    return entry.second;
}

// A reference into this thread's memo (map nodes do not move).
const CoreLoopPropertyEstimate& estimate_loop_property(const CoreMaterial& coreMaterial, LoopProperty property) {
    size_t generation;
    const auto& sourced = sourced_loop_members(property, generation);
    const auto materialType = coreMaterial.get_material();
    const auto composition = coreMaterial.get_material_composition();
    // The estimate is a function of the member list and of the material's (type, composition, loop family) alone;
    // the Roshen loop asks for it at every point of a magnetisation curve. Exact-key memo, per thread, per member-list
    // generation; only successful estimates are stored (the refusal below names the material, so it is rebuilt).
    using EstimateKey = std::tuple<size_t, LoopProperty, MaterialType, std::optional<MaterialComposition>, std::optional<NanocrystallineLoop>>;
    thread_local std::map<EstimateKey, CoreLoopPropertyEstimate> estimates;
    const auto loopFamilyKey = nanocrystalline_loop(coreMaterial);
    EstimateKey key{generation, property, materialType, composition, loopFamilyKey};
    if (auto found = estimates.find(key); found != estimates.end()) {
        return found->second;
    }
    const bool physicalComposition = composition && composition.value() != MaterialComposition::PROPRIETARY;
    const std::string typeName(magic_enum::enum_name(materialType));
    // Nanocrystalline: only members of the material's own loop family (no borrowing across families).
    const auto loopFamily = nanocrystalline_loop(coreMaterial);
    const std::string familyLabel = loop_family_label(loopFamily);

    CoreLoopPropertyEstimate estimate;
    std::optional<LoopPool> pool;
    std::string skipped;
    if (physicalComposition) {
        std::vector<const LoopMember*> members;
        for (const auto& member : sourced) {
            if (member.materialType == materialType && member.loopFamily == loopFamily && member.composition &&
                member.composition.value() == composition.value()) {
                members.push_back(&member);
            }
        }
        estimate.tier = CoreLoopPropertyTier::CLASS;
        estimate.pool = "class " + typeName + " " + std::string(magic_enum::enum_name(composition.value())) + familyLabel;
        pool = pool_median(members);
        if (!pool) {
            skipped = estimate.pool + " has no sourced member";
        }
    }
    else {
        skipped = composition ? "its composition is PROPRIETARY" : "it has no materialComposition";
    }
    if (!pool) {
        std::vector<const LoopMember*> members;
        for (const auto& member : sourced) {
            if (member.materialType == materialType && member.loopFamily == loopFamily) {
                members.push_back(&member);
            }
        }
        estimate.tier = CoreLoopPropertyTier::MATERIAL_TYPE;
        estimate.pool = "material type " + typeName + familyLabel;
        pool = pool_median(members);
    }
    if (!pool) {
        throw MaterialException(ErrorCode::MATERIAL_DATA_MISSING,
            "Core material '" + coreMaterial.get_name() + "' has no " + loop_property_name(property) + " in MAS, and no " +
            typeName + familyLabel + " material in MAS carries a sourced " + loop_property_name(property) +
            ": the loop of its class is unknown, so the Roshen hysteresis loop cannot be fitted.");
    }
    estimate.classValue = pool->value;
    estimate.sourceMaterials = pool->sourceMaterials;
    estimate.numberDocuments = pool->numberDocuments;
    estimate.numberDataPoints = pool->numberDataPoints;
    std::ostringstream method;
    method << "the median " << (property == LoopProperty::REMANENCE ? "Br/Bs " : "Hc ") << pool->value
           << (property == LoopProperty::REMANENCE ? " of its " : " A/m of its ") << pool->numberDataPoints
           << " sourced data points";
    if (property == LoopProperty::REMANENCE) {
        method << " times the material's own Bs";
    }
    if (!skipped.empty()) {
        method << " (" << skipped << ")";
    }
    estimate.method = method.str();
    return estimates.emplace(std::move(key), std::move(estimate)).first->second;
}

// One WARNING per material (manufacturer + name) and property, as the heat-conductivity estimate.
// The message is composed only for the warning that is emitted (it is the same text either way).
template <typename MessageBuilder>
void warn_loop_estimate_once(const CoreMaterial& coreMaterial, LoopProperty property, const MessageBuilder& message) {
    static std::mutex warnedMutex;
    static std::set<std::string> warned;
    std::lock_guard<std::mutex> lock(warnedMutex);
    if (warned.insert(coreMaterial.get_manufacturer_info().get_name() + "/" + coreMaterial.get_name() + "/" + loop_property_field(property)).second) {
        OM_WARNING_M("AmplitudePermeability", message());
    }
}

// The value resolve_loop_property returns, with its warning and its exceptions; estimateOut (if given) is pointed at
// the class estimate used, or nullptr for the material's own data. get_hysteresis_loop needs only the value.
double resolve_loop_property_value(const CoreMaterial& coreMaterial, double temperature, LoopProperty property,
                                   const CoreLoopPropertyEstimate** estimateOut = nullptr) {
    if (estimateOut) {
        *estimateOut = nullptr;
    }
    if (has_loop_property(coreMaterial, property)) {
        const double value = property == LoopProperty::REMANENCE ? Core::get_remanence(coreMaterial, temperature)
                                                                 : Core::get_coercive_force(coreMaterial, temperature);
        // A loop keeps less flux at H = 0 than it carries saturated. A material's own remanence at or above its
        // saturation flux density is a data error (MaterialValidator: MAT_REMANENCE_SATURATION), and the loop built
        // from it is not the material's: Fair-Rite 95 with 0.8 T at 100 C (for 0.08 T) gave a magnetisation curve
        // that topped out at 0.13 T at 60 C (Henry, 2026-10-24).
        if (property == LoopProperty::REMANENCE) {
            const double saturation = Core::get_magnetic_flux_density_saturation(coreMaterial, temperature, false);
            if (value >= saturation) {
                throw MaterialException(ErrorCode::MATERIAL_INVALID_PROPERTY,
                    "The data of material " + coreMaterial.get_name() + " is not physical: its remanence at " +
                    format_significant(temperature) + " C, " + format_significant(value) + " T, is not below its saturation flux "
                    "density there, " + format_significant(saturation) + " T, so no hysteresis loop can be built from it. Choose "
                    "another material, or correct this material's remanence");
            }
        }
        return value;
    }
    const auto& estimate = estimate_loop_property(coreMaterial, property);
    double value;
    if (property == LoopProperty::REMANENCE) {
        value = estimate.classValue * Core::get_magnetic_flux_density_saturation(coreMaterial, temperature, false);
    }
    else {
        value = Core::get_ferromagnetic_property_to_curie(coreMaterial, {{Defaults().ambientTemperature, estimate.classValue}},
                                                          temperature, "estimated coercive force");
    }
    warn_loop_estimate_once(coreMaterial, property, [&]() {
        std::ostringstream message;
        message << "Core material '" << coreMaterial.get_name() << "' has no " << loop_property_name(property)
                << " in MAS: ESTIMATED " << value << (property == LoopProperty::REMANENCE ? " T" : " A/m") << " at "
                << temperature << " C, tier " << magic_enum::enum_name(estimate.tier) << ", from the " << estimate.pool << ": "
                << estimate.sourceMaterials.size() << " sourced materials, " << estimate.numberDocuments << " documents, "
                << estimate.numberDataPoints << " independent data points, using " << estimate.method << ".";
        return message.str();
    });
    if (estimateOut) {
        *estimateOut = &estimate;
    }
    return value;
}

CoreLoopProperty resolve_loop_property(const CoreMaterial& coreMaterial, double temperature, LoopProperty property) {
    const CoreLoopPropertyEstimate* estimate = nullptr;
    const double value = resolve_loop_property_value(coreMaterial, temperature, property, &estimate);
    if (!estimate) {
        return CoreLoopProperty{value, CoreLoopPropertyTier::MATERIAL, std::nullopt};
    }
    return CoreLoopProperty{value, estimate->tier, *estimate};
}

} // namespace

CoreLoopProperty resolve_core_material_remanence(const CoreMaterial& coreMaterial, double temperature) {
    return resolve_loop_property(coreMaterial, temperature, LoopProperty::REMANENCE);
}

CoreLoopProperty resolve_core_material_coercive_force(const CoreMaterial& coreMaterial, double temperature) {
    return resolve_loop_property(coreMaterial, temperature, LoopProperty::COERCIVE_FORCE);
}

std::optional<double> AmplitudePermeability::get_amplitude_permeability(std::string coreMaterialName, std::optional<double> magneticFluxDensityPeak, std::optional<double> magneticFieldStrengthPeak, double temperature) {
    CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);
    return get_amplitude_permeability(coreMaterial, magneticFluxDensityPeak, magneticFieldStrengthPeak, temperature);
}

std::optional<double> AmplitudePermeability::get_amplitude_permeability(const CoreMaterial& coreMaterial, std::optional<double> magneticFluxDensityPeak, std::optional<double> magneticFieldStrengthPeak, double temperature) {
    // Calculate amplitude permeability using B-H loop slope
    // This uses the average slope of the loop, not the instantaneous derivative
    
    double amplitudePermeability = 1;
    if (!magneticFieldStrengthPeak && !magneticFluxDensityPeak) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Either H or B must be specified");
    }
    if (magneticFieldStrengthPeak) {
        auto [upperPath, lowerPath] = BHLoopRoshenModel().get_hysteresis_loop(coreMaterial, temperature, std::nullopt, magneticFieldStrengthPeak.value());
        if (upperPath.get_x_points().size() > 1){
            double dx = fabs(upperPath.get_x_points()[1] - upperPath.get_x_points()[0]);
            if (dx < 1e-12) {
                return std::nullopt;
            }
            amplitudePermeability = fabs(upperPath.get_y_points()[1] - upperPath.get_y_points()[0]) / dx / constants.vacuumPermeability;
        }
        else {
            return std::nullopt;
        }
    }
    else {
        auto [upperPath, lowerPath] = BHLoopRoshenModel().get_hysteresis_loop(coreMaterial, temperature, magneticFluxDensityPeak.value(), std::nullopt);
        if (upperPath.get_x_points().size() > 1){
            double dx = fabs(upperPath.get_x_points()[1] - upperPath.get_x_points()[0]);
            if (dx < 1e-12) {
                return std::nullopt;
            }
            amplitudePermeability = fabs(upperPath.get_y_points()[1] - upperPath.get_y_points()[0]) / dx / constants.vacuumPermeability;
        }
        else {
            return std::nullopt;
        }
    }
    // ABT #1951: a NaN or infinite slope is not an answer; callers test the optional, not the number.
    if (!std::isfinite(amplitudePermeability)) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
            "Amplitude permeability of " + coreMaterial.get_name() + " came out as " + std::to_string(amplitudePermeability) +
            " at " + (magneticFieldStrengthPeak ? "H peak " + std::to_string(magneticFieldStrengthPeak.value()) + " A/m"
                                                : "B peak " + std::to_string(magneticFluxDensityPeak.value()) + " T") +
            ", " + std::to_string(temperature) + " C");
    }
    return amplitudePermeability;
}

std::map<std::string, double> BHLoopRoshenModel::get_major_loop_parameters(double saturationMagneticFieldStrength,
                                                                           double saturationMagneticFluxDensity,
                                                                           double coerciveForce,
                                                                           double remanence) {
    std::map<std::string, double> majorLoopParameters;
    double a1 = 0;
    double b1 = 0;
    double b2 = 0;
    double Hc = coerciveForce;
    double H0 = saturationMagneticFieldStrength;
    double B0 = saturationMagneticFluxDensity;
    double H1 = 0;
    double B1 = remanence;
    double H2 = -saturationMagneticFieldStrength;
    double B2 = -saturationMagneticFluxDensity;
    b1 = (H0 / B0 + Hc / B0 - H1 / B1 - Hc / B1) / (H0 - H1);
    a1 = (Hc - B1 * b1 * Hc) / B1;
    b2 = (H2 + Hc - B2 * a1) / (B2 * fabs(H2 + Hc));
    majorLoopParameters["a1"] = a1;
    majorLoopParameters["b1"] = b1;
    majorLoopParameters["b2"] = b2;
    return majorLoopParameters;
}

std::pair<Curve2D, Curve2D> BHLoopRoshenModel::get_hysteresis_loop(std::string coreMaterialName, 
                                                                    double temperature, 
                                                                    std::optional<double> magneticFluxDensityPeak, 
                                                                    std::optional<double> magneticFieldStrengthPeak) {
    CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);
    return get_hysteresis_loop(coreMaterial, temperature, magneticFluxDensityPeak, magneticFieldStrengthPeak);
}

std::pair<Curve2D, Curve2D> BHLoopRoshenModel::get_hysteresis_loop(const CoreMaterial& coreMaterial, 
                                                                    double temperature, 
                                                                    std::optional<double> magneticFluxDensityPeak, 
                                                                    std::optional<double> magneticFieldStrengthPeak) {

    // ABT #1736: remanence and coercive force reach zero at the Curie temperature (see
    // Core::get_remanence); at and above it there is no loop to build.
    if (coreMaterial.get_curie_temperature() && temperature >= coreMaterial.get_curie_temperature().value()) {
        throw MaterialException(ErrorCode::MATERIAL_ABOVE_CURIE_TEMPERATURE,
            "Core material " + coreMaterial.get_name() + " has no hysteresis loop at " + std::to_string(temperature) +
            " C: it is not ferromagnetic at or above its Curie temperature of " +
            std::to_string(coreMaterial.get_curie_temperature().value()) + " C");
    }
    // ABT #1951: the Roshen loop is fitted through (Hsat, Bsat), (0, Br) and (-Hc, 0); its b1 holds
    // -Hc/Br. Core::get_remanence and Core::get_coercive_force answer NaN for an absent record
    // (MaterialValidator reads that), and that NaN went through the fit and out of
    // get_amplitude_permeability as a value: every Nanoperm grade at any H, and 523 of the 1073 MAS
    // materials. Absent is missing data, not zero. The resolvers give the material's own value, or
    // the flagged class estimate Alf approved (2026-10-18), or throw when the class is unknown.
    // Inside a magnetisation-curve build these four depend only on the build's material and temperature, and are
    // computed once per build (CurveBuildScope.h); elsewhere, every call.
    std::array<double, 4> loopInputs;
    auto* scopedInputs = curve_build::loop_inputs_for(coreMaterial, temperature);
    if (scopedInputs && scopedInputs->has_value()) {
        loopInputs = scopedInputs->value();
    }
    else {
        loopInputs[0] = Core::get_magnetic_field_strength_saturation(coreMaterial, temperature);
        loopInputs[1] = Core::get_magnetic_flux_density_saturation(coreMaterial, temperature, false);
        loopInputs[2] = resolve_loop_property_value(coreMaterial, temperature, LoopProperty::COERCIVE_FORCE);
        loopInputs[3] = resolve_loop_property_value(coreMaterial, temperature, LoopProperty::REMANENCE);
        if (scopedInputs) {
            *scopedInputs = loopInputs;
        }
    }
    double saturationMagneticFieldStrength = loopInputs[0];
    double saturationMagneticFluxDensity = loopInputs[1];
    double coerciveForce = loopInputs[2];
    double remanence = loopInputs[3];

    auto majorLoopParameters = get_major_loop_parameters(saturationMagneticFieldStrength, saturationMagneticFluxDensity,
                                                         coerciveForce, remanence);

    double a1 = majorLoopParameters["a1"];
    double b1 = majorLoopParameters["b1"];
    double b2 = majorLoopParameters["b2"];

    // OPTIMIZED: Use fewer points with adaptive spacing instead of fixed step
    std::vector<double> magneticFieldStrengthPoints;
    magneticFieldStrengthPoints.reserve(32);  // the critical, peak and bound points below

    // Critical points that must be included
    magneticFieldStrengthPoints.push_back(-saturationMagneticFieldStrength);
    magneticFieldStrengthPoints.push_back(-coerciveForce);
    magneticFieldStrengthPoints.push_back(0);
    magneticFieldStrengthPoints.push_back(coerciveForce);
    magneticFieldStrengthPoints.push_back(saturationMagneticFieldStrength);
    
    // Add points near origin for accurate permeability calculation
    // These are critical for slope calculation at low H
    double H_near_zero = saturationMagneticFieldStrength * 0.001;  // Very close to origin
    magneticFieldStrengthPoints.push_back(H_near_zero);
    magneticFieldStrengthPoints.push_back(-H_near_zero);
    magneticFieldStrengthPoints.push_back(coerciveForce * 0.1);
    magneticFieldStrengthPoints.push_back(-coerciveForce * 0.1);
    
    // Add points around the operating region if specified
    if (magneticFieldStrengthPeak) {
        double H_peak = magneticFieldStrengthPeak.value();
        // Add points near the peak for better resolution
        if (H_peak > 0) {
            magneticFieldStrengthPoints.push_back(H_peak * 0.001);  // Very close to 0
            magneticFieldStrengthPoints.push_back(H_peak * 0.01);   // Close to 0
            magneticFieldStrengthPoints.push_back(H_peak * 0.1);
            magneticFieldStrengthPoints.push_back(H_peak * 0.5);
            magneticFieldStrengthPoints.push_back(H_peak * 0.9);
            magneticFieldStrengthPoints.push_back(H_peak * 0.95);
            magneticFieldStrengthPoints.push_back(H_peak);
            magneticFieldStrengthPoints.push_back(-H_peak * 0.001);
            magneticFieldStrengthPoints.push_back(-H_peak * 0.01);
            magneticFieldStrengthPoints.push_back(-H_peak * 0.1);
            magneticFieldStrengthPoints.push_back(-H_peak * 0.5);
            magneticFieldStrengthPoints.push_back(-H_peak * 0.9);
            magneticFieldStrengthPoints.push_back(-H_peak * 0.95);
            magneticFieldStrengthPoints.push_back(-H_peak);
        }
    }
    
    // Sort and remove duplicates
    std::sort(magneticFieldStrengthPoints.begin(), magneticFieldStrengthPoints.end());
    auto last = std::unique(magneticFieldStrengthPoints.begin(), magneticFieldStrengthPoints.end(),
        [](double a, double b) { return fabs(a - b) < 1e-6; });
    magneticFieldStrengthPoints.erase(last, magneticFieldStrengthPoints.end());
    
    // Ensure bounds
    if (magneticFieldStrengthPoints.front() > -saturationMagneticFieldStrength) {
        magneticFieldStrengthPoints.insert(magneticFieldStrengthPoints.begin(), -saturationMagneticFieldStrength);
    }
    if (magneticFieldStrengthPoints.back() < saturationMagneticFieldStrength) {
        magneticFieldStrengthPoints.push_back(saturationMagneticFieldStrength);
    }

    auto bh_curve_half_loop = [&](double magneticFieldStrength, double a, double b) {
        return (magneticFieldStrength + coerciveForce) / (a + b * fabs(magneticFieldStrength + coerciveForce));
    };

    auto calculate_magnetic_flux_density = [&](double magneticFieldStrength, bool loop_is_upper = true) {
        double magneticFluxDensity;
        if (loop_is_upper) {
            if (magneticFieldStrength < -coerciveForce) {
                magneticFluxDensity = bh_curve_half_loop(magneticFieldStrength, a1, b2);
            }
            else {
                magneticFluxDensity = bh_curve_half_loop(magneticFieldStrength, a1, b1);
            }
        }
        else {
            if (magneticFieldStrength < coerciveForce) {
                magneticFluxDensity = -bh_curve_half_loop(-magneticFieldStrength, a1, b1);
            }
            else {
                magneticFluxDensity = -bh_curve_half_loop(-magneticFieldStrength, a1, b2);
            }
        }

        return magneticFluxDensity;
    };

    auto calculate_magnetic_flux_density_waveform = [&](const std::vector<double>& magneticFieldStrengthWaveform,
                                                  bool loop_is_upper = true) {
        std::vector<double> magneticFluxDensityWaveform;
        magneticFluxDensityWaveform.reserve(magneticFieldStrengthWaveform.size());
        for (auto& magneticFieldStrength : magneticFieldStrengthWaveform) {
            double magneticFluxDensity = calculate_magnetic_flux_density(magneticFieldStrength, loop_is_upper);
            magneticFluxDensityWaveform.push_back(magneticFluxDensity);
        }

        return magneticFluxDensityWaveform;
    };

    std::vector<double> upperMagneticFluxDensityWaveform =
        calculate_magnetic_flux_density_waveform(magneticFieldStrengthPoints, true);
    std::vector<double> lowerMagneticFluxDensityWaveform =
        calculate_magnetic_flux_density_waveform(magneticFieldStrengthPoints, false);

    Curve2D upperPath;
    Curve2D lowerPath;
    upperPath.get_mutable_x_points().reserve(magneticFieldStrengthPoints.size());
    upperPath.get_mutable_y_points().reserve(magneticFieldStrengthPoints.size());
    lowerPath.get_mutable_x_points().reserve(magneticFieldStrengthPoints.size());
    lowerPath.get_mutable_y_points().reserve(magneticFieldStrengthPoints.size());

    if (magneticFieldStrengthPeak) {
        // OPTIMIZED: Find index with binary search instead of min_element
        double H_peak = magneticFieldStrengthPeak.value();
        size_t indexOfDesiredHPeak = 0;
        double minDiff = fabs(magneticFieldStrengthPoints[0] - H_peak);
        for (size_t i = 1; i < magneticFieldStrengthPoints.size(); ++i) {
            double diff = fabs(magneticFieldStrengthPoints[i] - H_peak);
            if (diff < minDiff) {
                minDiff = diff;
                indexOfDesiredHPeak = i;
            }
        }
        
        double differenceAtDesiredHPeak = fabs(upperMagneticFluxDensityWaveform[indexOfDesiredHPeak] - 
                                                lowerMagneticFluxDensityWaveform[indexOfDesiredHPeak]);

        for (size_t i = 0; i < upperMagneticFluxDensityWaveform.size(); i++) {
            upperMagneticFluxDensityWaveform[i] -= differenceAtDesiredHPeak / 2;
        }
        for (size_t i = 0; i < lowerMagneticFluxDensityWaveform.size(); i++) {
            lowerMagneticFluxDensityWaveform[i] += differenceAtDesiredHPeak / 2;
        }

        for (size_t i = 0; i < upperMagneticFluxDensityWaveform.size(); i++) {
            auto elemB = upperMagneticFluxDensityWaveform[i];
            auto elemH = magneticFieldStrengthPoints[i];
            if (elemH <= magneticFieldStrengthPeak.value() && elemH >= -magneticFieldStrengthPeak.value()) {
                upperPath.get_mutable_x_points().push_back(elemH);
                upperPath.get_mutable_y_points().push_back(elemB);
            }
        }

        for (size_t i = 0; i < lowerMagneticFluxDensityWaveform.size(); i++) {
            auto elemB = lowerMagneticFluxDensityWaveform[i];
            auto elemH = magneticFieldStrengthPoints[i];
            if (elemH <= magneticFieldStrengthPeak.value() && elemH >= -magneticFieldStrengthPeak.value()) {
                lowerPath.get_mutable_x_points().push_back(elemH);
                lowerPath.get_mutable_y_points().push_back(elemB);
            }
        }
    }
    else if (magneticFluxDensityPeak) {
        double magneticFluxDensityDifference = magneticFluxDensityPeak.value();
        size_t timeout = 0;
        double abs_tol = 0.001;
        
        // OPTIMIZED: Pre-allocate difference vector
        std::vector<double> difference;
        difference.reserve(upperMagneticFluxDensityWaveform.size());
        
        while (fabs(magneticFluxDensityDifference) > abs_tol && timeout < 10) {
            difference.clear();
            for (size_t i = 0; i < upperMagneticFluxDensityWaveform.size(); i++) {
                difference.push_back(fabs(upperMagneticFluxDensityWaveform[i] - lowerMagneticFluxDensityWaveform[i]));
            }
            
            size_t mininumPosition =
                std::distance(std::begin(difference), std::min_element(std::begin(difference), std::end(difference)));
            magneticFluxDensityDifference =
                fabs(upperMagneticFluxDensityWaveform[mininumPosition]) - magneticFluxDensityPeak.value();

            for (size_t i = 0; i < upperMagneticFluxDensityWaveform.size(); i++) {
                upperMagneticFluxDensityWaveform[i] -= magneticFluxDensityDifference / 16;
            }
            for (size_t i = 0; i < lowerMagneticFluxDensityWaveform.size(); i++) {
                lowerMagneticFluxDensityWaveform[i] += magneticFluxDensityDifference / 16;
            }

            timeout++;
            abs_tol += timeout * 0.0001;
        }

        for (size_t i = 0; i < upperMagneticFluxDensityWaveform.size(); i++) {
            auto elemB = upperMagneticFluxDensityWaveform[i];
            auto elemH = magneticFieldStrengthPoints[i];
            if (elemB <= magneticFluxDensityPeak.value() && elemB >= -magneticFluxDensityPeak.value()) {
                upperPath.get_mutable_x_points().push_back(elemH);
                upperPath.get_mutable_y_points().push_back(elemB);
            }
        }

        for (size_t i = 0; i < lowerMagneticFluxDensityWaveform.size(); i++) {
            auto elemB = lowerMagneticFluxDensityWaveform[i];
            auto elemH = magneticFieldStrengthPoints[i];
            if (elemB <= magneticFluxDensityPeak.value() && elemB >= -magneticFluxDensityPeak.value()) {
                lowerPath.get_mutable_x_points().push_back(elemH);
                lowerPath.get_mutable_y_points().push_back(elemB);
            }
        }
    }

    return {upperPath, lowerPath};
}

} // namespace OpenMagnetics
