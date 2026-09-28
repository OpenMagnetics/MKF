// SPDX-License-Identifier: MIT
// "Magnetic Blade Runner": the physics validator for MAS core materials.
//
// Given one MAS core-material record (one line of core_materials.ndjson, as JSON), decide
// whether its data is physically believable and say why. The finding format is Blade Runner's
// (PSMA/TAS/validator, tas::Finding / tas::Verdict): a stable rule code, a severity, the record
// it belongs to, a message naming the offending quantity, the value and the bound it broke.
//
//   IMPOSSIBLE  the data contradicts physics (a record with one is INVALID; fails MAS CI)
//   SUSPICIOUS  implausible for the material class, or not evaluable where it will be used
//   WARNING     a data-policy gap that says nothing about the physics (MAS RFC 0011 provenance)
//
// Every number the validator judges is produced by MKF's own models (CoreLosses, Core,
// InitialPermeability): nothing here re-implements a loss formula. The class envelopes are
// PUBLISHED datasheet numbers, each cited where it is defined (MaterialValidatorEnvelopes.cpp).
//
// A check that cannot run on a record says so in `skipped`, with the reason. A check that
// examines nothing looks exactly like a check that found nothing wrong (TAS ABT #387), so a
// material missing the data a check needs is reported, never silently passed and never given
// a default.
#pragma once

#include <nlohmann/json.hpp>

#include <map>
#include <optional>
#include <string>
#include <vector>

namespace OpenMagnetics {

using json = nlohmann::json;

enum class MaterialFindingSeverity { WARNING, SUSPICIOUS, IMPOSSIBLE };

const char* to_string(MaterialFindingSeverity severity);

struct MaterialFinding {
    std::string code;          // stable rule id, e.g. "MAT_LOSS_ENVELOPE"
    MaterialFindingSeverity severity = MaterialFindingSeverity::SUSPICIOUS;
    std::string component = "coreMaterial";
    std::string reference;     // the material name
    std::string message;       // names the offending quantity and the operating point
    double value = 0.0;        // the computed quantity that tripped the rule
    double threshold = 0.0;    // the bound it violated (0 if not applicable)
};

struct MaterialVerdict {
    std::string reference;                 // the material name
    bool valid = true;                     // false iff any finding is IMPOSSIBLE
    std::vector<MaterialFinding> findings; // every rule that fired
    std::vector<std::string> skipped;      // "<RULE>: <reason>" for every check that could not run
    std::string materialClass;             // the class the envelopes were taken from
};

// One sinusoidal volumetric-loss measurement point, read from a record's points array
// (core_materials.ndjson or advanced_core_materials.ndjson).
struct MaterialLossPoint {
    double frequency;                // Hz
    double magneticFluxDensityPeak;  // T
    double temperature;              // C
    double volumetricLosses;         // W/m^3
    std::string origin;
};

// A reference operating point of a material class, with the best and worst PUBLISHED grades there.
struct MaterialLossAnchor {
    std::string manufacturer;
    std::string grade;
    double volumetricLosses;  // W/m^3, sinusoidal
    std::string source;       // document, page/table/figure
};

struct MaterialLossReferencePoint {
    double frequency;                // Hz
    double magneticFluxDensityPeak;  // T
    double temperature;              // C
    // The point that defines the class (e.g. 100 kHz for MnZn power ferrite). A loss model whose
    // fitted span does not reach it is MAT_LOSS_OUT_OF_RANGE: a consumer using the grade at the
    // class's standard operating point would be extrapolating. Non-primary points are judged when
    // they lie inside a declared fitted span (reported in `skipped` otherwise), and on the LOW side
    // only: their worst anchor covers only the grades marketed at that point.
    bool primary;
    MaterialLossAnchor best;
    MaterialLossAnchor worst;
};

struct MaterialClassBounds {
    std::string name;                  // e.g. "MnZn power ferrite"
    std::vector<MaterialLossReferencePoint> lossReferencePoints;  // empty: no published anchor
    std::string lossAnchorGap;         // why there is no loss envelope, when there is none
    // WEAK class: the anchors are thin (curve readings, fits, one maker, unstated temperature),
    // so a breach is SUSPICIOUS only, never IMPOSSIBLE.
    bool weak = false;
    std::string anchorQuality;         // one line on how good the anchors are
    // Plausibility of the static properties (SUSPICIOUS outside), each with its source.
    std::optional<double> minimumInitialPermeability;
    std::optional<double> maximumInitialPermeability;
    std::string initialPermeabilitySource;
    std::optional<double> maximumSaturationFluxDensity;  // T at 25 C
    std::string saturationSource;
    std::optional<double> minimumCurieTemperature;  // C
    std::optional<double> maximumCurieTemperature;  // C
    std::string curieSource;
};

// The class a record belongs to, decided ONLY from its MAS type / materialComposition /
// application fields; nullptr when those fields do not determine a class (never guessed).
const MaterialClassBounds* classify_material(const json& materialRecord);

class MaterialValidator {
  public:
    // Validate one core-material record. `extraLossPoints` are measurement points kept outside
    // the record (advanced_core_materials.ndjson); the record's own points are read from it.
    // Never throws on bad DATA: a record MKF cannot read is itself a finding (MAT_PARSE).
    MaterialVerdict validate(const json& materialRecord,
                             const std::vector<MaterialLossPoint>& extraLossPoints = {}) const;

    // Validate every record of a core_materials.ndjson file, merging the loss points of an
    // advanced_core_materials.ndjson file when given (the advanced file is streamed; only its
    // sinusoidal loss points are kept).
    std::vector<MaterialVerdict> validate_catalogue(const std::string& coreMaterialsPath,
                                                    const std::optional<std::string>& advancedCoreMaterialsPath = std::nullopt) const;
    // Same, over the core_materials.ndjson MKF was built with (embedded), no advanced points.
    std::vector<MaterialVerdict> validate_embedded_catalogue() const;
    std::vector<MaterialVerdict> validate_catalogue_text(const std::string& coreMaterialsNdjson,
                                                         const std::map<std::string, std::vector<MaterialLossPoint>>& extraLossPoints = {}) const;

    // Sinusoidal loss points of a record's volumetricLosses points arrays (any shape-family key).
    static std::vector<MaterialLossPoint> read_loss_points(const json& materialRecord);
    static std::map<std::string, std::vector<MaterialLossPoint>> read_advanced_loss_points(const std::string& advancedCoreMaterialsPath);

    // All rule codes this validator can emit.
    static std::vector<std::string> check_codes();
    // The class table (envelopes and plausibility bounds, with their sources).
    static const std::vector<MaterialClassBounds>& class_table();
    static json class_table_json();
};

void to_json(json& j, const MaterialFinding& finding);
void to_json(json& j, const MaterialVerdict& verdict);

// Per-rule / per-severity counts and the IMPOSSIBLE + SUSPICIOUS list of a sweep.
json summarize_material_verdicts(const std::vector<MaterialVerdict>& verdicts);

} // namespace OpenMagnetics
