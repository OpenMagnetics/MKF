#pragma once
#include "constructive_models/Wire.h"
#include "constructive_models/Core.h"
#include "constructive_models/Bobbin.h"
#include "constructive_models/InsulationMaterial.h"
#include <MAS.hpp>
#include "json.hpp"

#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace OpenMagnetics {

using nlohmann::json;

// Per-call container for a user-supplied catalog of cores, shapes, materials,
// wires, bobbins, insulation materials, and wire materials. The caller
// constructs a LibraryContext, fills it via loadFrom*(...), then passes a
// pointer to an adviser's public entry point. While the adviser runs, a
// RAII guard temporarily replaces the process-global *Database maps with
// the per-kind data from this context (Merge: union with built-ins, custom
// wins on name clash; Replace: only kinds the user provided are wiped, the
// rest fall through to built-ins).
//
// Orchestrating-thread only: Scope swaps the SHARED global catalogs, so it
// must be applied before set_databases_frozen(true) and destroyed after
// unfreezing. Constructing a Scope while databases are frozen throws
// (ABT #113 — see the THREAD-SAFETY CONTRACT in support/Utils.h).
class LibraryContext {
public:
    enum class LoadMode { Merge, Replace };

    LibraryContext() = default;

    LoadMode mode() const { return _mode; }
    void setMode(LoadMode m) { _mode = m; }

    // Parse a JSON document with optional top-level keys: "cores",
    // "coreMaterials", "coreShapes", "wires", "bobbins",
    // "insulationMaterials", "wireMaterials". Each is either a map keyed by
    // name or an array of objects whose "name" field is the key. Missing
    // keys are simply not populated (caller can still provide them via
    // another call). Throws with a kind-qualified message on parse error;
    // never silently substitutes a default.
    void loadFromJson(const json& data, LoadMode mode);
    void loadFromString(const std::string& jsonText, LoadMode mode);
    void loadFromFile(const std::string& path, LoadMode mode);

    void clear();
    bool empty() const;

    // RAII helper: on construction, swaps the global *Database maps with
    // a view derived from this context (Merge or Replace per `_mode`); on
    // destruction, restores them. Use via LibraryContext::applyScoped().
    class Scope {
    public:
        Scope() = default;
        Scope(const LibraryContext* ctx);
        ~Scope();
        Scope(const Scope&) = delete;
        Scope& operator=(const Scope&) = delete;
        Scope(Scope&& other) noexcept;
        Scope& operator=(Scope&& other) noexcept;
    private:
        bool _active = false;
        std::vector<Core> _coreBackup;
        std::map<std::string, MAS::CoreMaterial> _coreMaterialBackup;
        std::map<std::string, MAS::CoreShape> _coreShapeBackup;
        std::vector<MAS::CoreShapeFamily> _coreShapeFamilyBackup;
        std::map<std::string, Wire> _wireBackup;
        std::map<std::string, Bobbin> _bobbinBackup;
        std::map<std::string, InsulationMaterial> _insulationMaterialBackup;
        std::map<std::string, MAS::WireMaterial> _wireMaterialBackup;
    };

    Scope applyScoped() const { return Scope(this); }

    // Accessors (read-only) used by tests/bindings.
    const std::optional<std::vector<Core>>& cores() const { return _cores; }
    const std::optional<std::map<std::string, MAS::CoreMaterial>>& coreMaterials() const { return _coreMaterials; }
    const std::optional<std::map<std::string, MAS::CoreShape>>& coreShapes() const { return _coreShapes; }
    const std::optional<std::map<std::string, Wire>>& wires() const { return _wires; }
    const std::optional<std::map<std::string, Bobbin>>& bobbins() const { return _bobbins; }
    const std::optional<std::map<std::string, InsulationMaterial>>& insulationMaterials() const { return _insulationMaterials; }
    const std::optional<std::map<std::string, MAS::WireMaterial>>& wireMaterials() const { return _wireMaterials; }

private:
    LoadMode _mode = LoadMode::Merge;
    std::optional<std::vector<Core>>                                _cores;
    std::optional<std::map<std::string, MAS::CoreMaterial>>         _coreMaterials;
    std::optional<std::map<std::string, MAS::CoreShape>>            _coreShapes;
    std::optional<std::map<std::string, Wire>>                      _wires;
    std::optional<std::map<std::string, Bobbin>>                    _bobbins;
    std::optional<std::map<std::string, InsulationMaterial>>        _insulationMaterials;
    std::optional<std::map<std::string, MAS::WireMaterial>>         _wireMaterials;
};

// Convenience: type filter (allow + block lists). Empty allow = no
// allowlist constraint. Empty block = no blocklist constraint. A value
// passes iff (allow.empty() || allow.contains(v)) && !block.contains(v).
struct TypeFilterSet {
    std::set<std::string> allowed;
    std::set<std::string> blocked;

    bool empty() const { return allowed.empty() && blocked.empty(); }
    bool accepts(const std::string& v) const {
        if (!allowed.empty() && allowed.find(v) == allowed.end()) return false;
        if (blocked.find(v) != blocked.end()) return false;
        return true;
    }
};

// Bag of type-based hard constraints applied as a pre-filter at the top
// of each adviser pipeline. Values are MAS enum names as strings
// (case-insensitive, see acceptsCoreShapeFamily / acceptsMaterialType /
// acceptsWireType for the exact matching rules).
struct AdviserConstraints {
    TypeFilterSet shapeFamily;
    TypeFilterSet coreMaterialType;
    TypeFilterSet wireType;

    bool empty() const {
        return shapeFamily.empty() && coreMaterialType.empty() && wireType.empty();
    }
};

// Predicates: case-insensitive match against the magic_enum name of the
// MAS enum value (e.g. CoreShapeFamily::ETD matches "ETD"/"etd"/"Etd").
bool acceptsCoreShapeFamily(const TypeFilterSet& f, MAS::CoreShapeFamily family);
bool acceptsCoreMaterialType(const TypeFilterSet& f, MAS::MaterialType type);
bool acceptsWireType(const TypeFilterSet& f, MAS::WireType type);

// Apply the constraints to a candidate core list. Removes cores whose
// shape family / material type fails the filter. Returns a new vector.
std::vector<Core> filterCoresByConstraints(const std::vector<Core>& cores,
                                            const AdviserConstraints& constraints);

// Apply the wireType constraint to a candidate wire list. Returns a new
// map containing only accepted wires.
std::map<std::string, Wire> filterWiresByConstraints(
    const std::map<std::string, Wire>& wires,
    const AdviserConstraints& constraints);

} // namespace OpenMagnetics
