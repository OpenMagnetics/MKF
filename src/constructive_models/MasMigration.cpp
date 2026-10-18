// Backwards-compatibility helper for loading pre-1.0 MAS documents.
//
// See MasMigration.h for the contract. The MAPPING table below is the
// C++ port of MAS/scripts/migrate-to-1.0.py.

#include "constructive_models/MasMigration.h"
#include "Defaults.h"

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMagnetics::compat {

namespace {

// An old spelling is rewritten only when it is the value of one of the listed keys (or an element
// of an array under one): the keys whose CURRENT schema enum carries the new spelling. Matching by
// value alone rewrote free text too, e.g. a winding named "Forward" or a bobbin named "Basic"
// (ABT #1625; the same bug as MAS #1392 in migrate-to-1.0.py). The key lists were derived from the
// MAS schemas' properties whose enum contains the new value; subApplication (now a PEAS string)
// and the three topologies MAS later dropped keep the key they had before 1.0.
struct EnumMigration {
    std::string newValue;
    std::vector<std::string> keys;
};

const std::unordered_map<std::string, EnumMigration>& mapping() {
    static const std::unordered_map<std::string, EnumMigration> m = {
        // insulationType
        {"Functional", {"functional", {"insulationGrade", "insulationType", "type"}}},
        {"Basic", {"basic", {"insulationGrade", "insulationType"}}},
        {"Supplementary", {"supplementary", {"insulationGrade", "insulationType"}}},
        {"Double", {"double", {"insulationGrade", "insulationType"}}},
        {"Reinforced", {"reinforced", {"insulationGrade", "insulationType"}}},
        // wiringTechnology
        {"Wound", {"wound", {"type", "wiringTechnology"}}},
        {"Printed", {"printed", {"type", "wiringTechnology"}}},
        {"Stamped", {"stamped", {"type", "wiringTechnology"}}},
        {"Deposition", {"deposition", {"type", "wiringTechnology"}}},
        // application
        {"Power", {"power", {"application", "function"}}},
        {"Signal Processing", {"signalProcessing", {"application"}}},
        {"Interference Suppression", {"interferenceSuppression", {"application"}}},
        // subApplication
        {"Power Filtering", {"powerFiltering", {"subApplication"}}},
        {"Transforming", {"transforming", {"subApplication"}}},
        {"Isolation", {"isolation", {"subApplication"}}},
        {"Common Mode Noise Filtering", {"commonModeNoiseFiltering", {"subApplication"}}},
        {"Differential Mode Noise Filtering", {"differentialModeNoiseFiltering", {"subApplication"}}},
        // connectionType
        {"Pin", {"pin", {"assemblyType", "mounting", "terminalType", "type"}}},
        {"Screw", {"screw", {"assemblyType", "mounting", "terminalType", "type"}}},
        {"SMT", {"smt", {"assemblyType", "mounting", "terminalType", "type"}}},
        {"Flying Lead", {"flyingLead", {"assemblyType", "mounting", "terminalType", "type"}}},
        {"THT", {"tht", {"assemblyType", "mounting", "terminalType", "type"}}},
        {"PCB Pad", {"pcbPad", {"assemblyType", "mounting", "terminalType", "type"}}},
        // coilAlignment
        {"inner or top", {"innerOrTop", {"layersAlignment", "sectionsAlignment", "turnsAlignment"}}},
        {"outer or bottom", {"outerOrBottom", {"layersAlignment", "sectionsAlignment", "turnsAlignment"}}},
        // coreType
        {"two-piece set", {"twoPieceSet", {"type"}}},
        {"piece and plate", {"pieceAndPlate", {"type"}}},
        {"closed shape", {"closedShape", {"type"}}},
        // piece type
        {"half set", {"halfSet", {"type"}}},
        // DMC configuration
        {"SINGLE_PHASE", {"singlePhase", {"configuration"}}},
        {"THREE_PHASE", {"threePhase", {"configuration"}}},
        {"THREE_PHASE_WITH_NEUTRAL", {"threePhaseWithNeutral", {"configuration"}}},
        // waveformLabel
        {"Custom", {"custom", {"ancillaryLabel", "label", "pattern", "type", "waveformLabel"}}},
        {"Triangular", {"triangular", {"ancillaryLabel", "label", "waveformLabel"}}},
        {"Sinusoidal", {"sinusoidal", {"ancillaryLabel", "label", "waveformLabel"}}},
        {"Rectangular", {"rectangular", {"ancillaryLabel", "columnShape", "crossSectionalShape", "label", "shape", "type", "waveformLabel"}}},
        {"Unipolar Rectangular", {"unipolarRectangular", {"ancillaryLabel", "label", "waveformLabel"}}},
        {"Unipolar Triangular", {"unipolarTriangular", {"ancillaryLabel", "label", "waveformLabel"}}},
        {"Bipolar Rectangular", {"bipolarRectangular", {"ancillaryLabel", "label", "waveformLabel"}}},
        {"Bipolar Triangular", {"bipolarTriangular", {"ancillaryLabel", "label", "waveformLabel"}}},
        {"Flyback Primary", {"flybackPrimary", {"ancillaryLabel", "label", "waveformLabel"}}},
        {"Flyback Secondary", {"flybackSecondary", {"ancillaryLabel", "label", "waveformLabel"}}},
        {"Rectangular With Deadtime", {"rectangularWithDeadtime", {"ancillaryLabel", "label", "waveformLabel"}}},
        {"Flyback Secondary With Deadtime", {"flybackSecondaryWithDeadtime", {"ancillaryLabel", "label", "waveformLabel"}}},
        {"Rectangular DCM", {"rectangularDCM", {"ancillaryLabel", "label", "waveformLabel"}}},
        {"Secondary Rectangular", {"secondaryRectangular", {"ancillaryLabel", "label", "waveformLabel"}}},
        {"Secondary Rectangular With Deadtime", {"secondaryRectangularWithDeadtime", {"ancillaryLabel", "label", "waveformLabel"}}},
        {"Triangular With Deadtime", {"triangularWithDeadtime", {"ancillaryLabel", "label", "waveformLabel"}}},
        // topology
        {"Buck Converter", {"buckConverter", {"topology"}}},
        {"Boost Converter", {"boostConverter", {"topology"}}},
        {"Inverting Buck-Boost Converter", {"invertingBuckBoostConverter", {"topology"}}},
        {"Cuk Converter", {"cukConverter", {"topology"}}},
        {"SEPIC Converter", {"sepicConverter", {"topology"}}},
        {"Zeta Converter", {"zetaConverter", {"topology"}}},
        {"Weinberg Converter", {"weinbergConverter", {"topology"}}},
        {"Flyback Converter", {"flybackConverter", {"topology"}}},
        {"Active Clamp Forward Converter", {"activeClampForwardConverter", {"topology"}}},
        {"Single Switch Forward Converter", {"singleSwitchForwardConverter", {"topology"}}},
        {"Two Switch Forward Converter", {"twoSwitchForwardConverter", {"topology"}}},
        {"Push-Pull Converter", {"pushPullConverter", {"topology"}}},
        {"Weinberg Converter", {"weinbergConverter", {"topology"}}},
        {"4-Switch Buck-Boost Converter", {"fourSwitchBuckBoostConverter", {"topology"}}},
        {"Half-Bridge Converter", {"halfBridgeConverter", {"topology"}}},
        {"Full-Bridge Converter", {"fullBridgeConverter", {"topology"}}},
        {"Phase-Shifted Full-Bridge Converter", {"phaseShiftedFullBridgeConverter", {"topology"}}},
        {"Phase-Shifted Half-Bridge Converter", {"phaseShiftedHalfBridgeConverter", {"topology"}}},
        {"Current Transformer", {"currentTransformer", {"topology"}}},
        {"Isolated Buck Converter", {"isolatedBuckConverter", {"topology"}}},
        {"Isolated Buck-Boost Converter", {"isolatedBuckBoostConverter", {"topology"}}},
        {"Dual Active Bridge Converter", {"dualActiveBridgeConverter", {"topology"}}},
        {"LLC Resonant Converter", {"llcResonantConverter", {"topology"}}},
        {"CLLC Resonant Converter", {"cllcResonantConverter", {"topology"}}},
        {"CLLLC Resonant Converter", {"clllcResonantConverter", {"topology"}}},
        {"Common Mode Choke", {"commonModeChoke", {"subtype", "topology"}}},
        {"Differential Mode Choke", {"differentialModeChoke", {"topology"}}},
        {"Power Factor Correction", {"powerFactorCorrection", {"topology"}}},
        // flybackModes
        {"Continuous Conduction Mode", {"continuousConductionMode", {"mode"}}},
        {"Discontinuous Conduction Mode", {"discontinuousConductionMode", {"mode"}}},
        {"Quasi Resonant Mode", {"quasiResonantMode", {"mode"}}},
        {"Boundary Mode Operation", {"boundaryModeOperation", {"mode"}}},
        // pfcModes
        {"Critical Conduction Mode", {"criticalConductionMode", {"mode"}}},
        {"Transition Mode", {"transitionMode", {"mode"}}},
        // bridge types
        {"Half Bridge", {"halfBridge", {"bridgeType", "bridgeTypePrimary", "bridgeTypeSecondary"}}},
        {"Full Bridge", {"fullBridge", {"bridgeType", "bridgeTypePrimary", "bridgeTypeSecondary", "rectifierType"}}},
        {"Center Tapped", {"centerTapped", {"rectifierType"}}},
        {"Current Doubler", {"currentDoubler", {"rectifierType"}}},
        // cllcPowerFlow
        {"Forward", {"forward", {"powerFlow", "powerFlowDirection"}}},
        {"Reverse", {"reverse", {"powerFlow", "powerFlowDirection"}}},
        // market
        {"Medical", {"medical", {"market"}}},
        {"Commercial", {"commercial", {"market", "type"}}},
        {"Industrial", {"industrial", {"market"}}},
        {"Military", {"military", {"market"}}},
        {"Space", {"space", {"market"}}},
        // materialComposition
        {"Carbonyl Iron", {"carbonylIron", {"materialComposition"}}},
        {"Iron", {"iron", {"materialComposition"}}},
        {"Proprietary", {"proprietary", {"materialComposition"}}},
        // cti
        {"Group I", {"groupI", {"cti"}}},
        {"Group II", {"groupII", {"cti"}}},
        {"Group IIIA", {"groupIIIA", {"cti"}}},
        {"Group IIIB", {"groupIIIB", {"cti"}}},
        // pollutionDegree (RFC 0008)
        {"P1", {"PD1", {"pollutionDegree"}}},
        {"P2", {"PD2", {"pollutionDegree"}}},
        {"P3", {"PD3", {"pollutionDegree"}}},
        // overvoltageCategory (RFC 0008)
        {"OVC-I", {"I", {"overvoltageCategory"}}},
        {"OVC-II", {"II", {"overvoltageCategory"}}},
        {"OVC-III", {"III", {"overvoltageCategory"}}},
        {"OVC-IV", {"IV", {"overvoltageCategory"}}},
        // shape family
        {"planar e", {"planarE", {"family"}}},
        {"planar er", {"planarER", {"family"}}},
        {"planar el", {"planarEL", {"family"}}},
    };
    return m;
}

// parentKey: the object key this value sits under (array elements inherit their array's key);
// nullptr at the root of the walk.
void migrate(nlohmann::json& j, const std::string* parentKey) {
    if (j.is_string()) {
        if (parentKey == nullptr) {
            return;
        }
        const auto& m = mapping();
        auto it = m.find(j.get_ref<const std::string&>());
        if (it != m.end() &&
            std::find(it->second.keys.begin(), it->second.keys.end(), *parentKey) != it->second.keys.end()) {
            j = it->second.newValue;
        }
    } else if (j.is_object()) {
        // Core spacers written by older MKF, both rewritten to the dielectric MKF now stamps --
        // Defaults().defaultSpacerMaterial, PET, Alf's ruling (2026-09-13, ABT #1200):
        //   1. before ABT #1170: {"type": "spacer", "material": "plastic"} with no
        //      insulationMaterial -- the piece-only `material` key, which core/spacer.json forbids
        //      (additionalProperties: false), and no insulationMaterial, which it requires;
        //   2. between ABT #1170 and #1200: {"type": "spacer", "insulationMaterial": "plastic"} --
        //      schema-valid, but "plastic" is not a name in data/insulation_materials.ndjson. The
        //      same value can arrive as an inline record, {"insulationMaterial": {"name": "plastic",
        //      ...}}, once a design has been round-tripped through a serialiser; that too is mapped,
        //      but ONLY when the record carries no relativePermittivity of its own. A record that does
        //      is real data -- the gap-dielectric resolver uses an inline record's own values -- and
        //      overriding it with PET would be exactly the guess the no-fallbacks rule forbids.
        // Matched on those EXACT legacy forms and nothing broader: any other spacer material, or a
        // spacer missing one, is a malformed record, is left untouched, and throws when its
        // dielectric is needed.
        const auto type = j.find("type");
        if (type != j.end() && type->is_string() && type->get<std::string>() == "spacer") {
            const auto material = j.find("material");
            const auto insulationMaterial = j.find("insulationMaterial");
            const bool legacyPieceMaterial = material != j.end() && material->is_string() &&
                                             material->get<std::string>() == "plastic" &&
                                             (insulationMaterial == j.end() || insulationMaterial->is_null());
            const bool legacyInsulationMaterial = insulationMaterial != j.end() && insulationMaterial->is_string() &&
                                                  insulationMaterial->get<std::string>() == "plastic";
            const bool legacyInsulationRecord = insulationMaterial != j.end() && insulationMaterial->is_object() &&
                                                insulationMaterial->contains("name") &&
                                                insulationMaterial->at("name").is_string() &&
                                                insulationMaterial->at("name").get<std::string>() == "plastic" &&
                                                (!insulationMaterial->contains("relativePermittivity") ||
                                                 insulationMaterial->at("relativePermittivity").is_null());
            if (legacyPieceMaterial) {
                j["insulationMaterial"] = Defaults().defaultSpacerMaterial;
                j.erase("material");
            }
            else if (legacyInsulationMaterial || legacyInsulationRecord) {
                j["insulationMaterial"] = Defaults().defaultSpacerMaterial;
            }
        }
        for (auto it = j.begin(); it != j.end(); ++it) {
            // Current schema: distributorInfo.cost is a {value, currency} object (CurrencyAmount).
            // Normalize legacy forms UP to it so old data still loads:
            //   bare number   -> {value: <number>, currency: "USD"}  (pre-CurrencyAmount / PEAS-number era)
            //   string        -> null
            //   object / null -> pass through UNCHANGED
            // The object form is the norm, so the common path does only a couple of cheap type checks
            // and never copies/rewrites — the fast path is not slowed. (Previously this FLATTENED the
            // {value,currency} object back to a number, which then crashed from_json(CurrencyAmount).)
            if (it.key() == "cost") {
                if (it.value().is_number()) {
                    const double costValue = it.value().get<double>();
                    it.value() = nlohmann::json{{"value", costValue}, {"currency", "USD"}};
                } else if (it.value().is_string()) {
                    it.value() = nlohmann::json(nullptr);
                }
                // object and null forms pass through unchanged
            } else {
                migrate(it.value(), &it.key());
            }
        }
    } else if (j.is_array()) {
        for (auto& e : j) {
            migrate(e, parentKey);
        }
    }
}

}  // namespace

void migrate_pre_1_0(nlohmann::json& j) {
    migrate(j, nullptr);
}

}  // namespace OpenMagnetics::compat
