#include "support/LibraryContext.h"
#include "support/Utils.h"
#include "support/Exceptions.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <magic_enum.hpp>
#include <sstream>
#include <stdexcept>

namespace OpenMagnetics {

namespace {

// Iterate a JSON node that may be either an object map (name -> entry) or
// an array of entries. Calls cb(name, entryJson).
template <typename Fn>
void forEachEntry(const json& node, const std::string& kind, Fn&& cb) {
    if (node.is_null()) return;
    if (node.is_array()) {
        for (auto& el : node) {
            if (!el.contains("name") || !el["name"].is_string()) {
                throw std::runtime_error(
                    "LibraryContext: entry in '" + kind + "' has no string 'name' field");
            }
            cb(el["name"].get<std::string>(), el);
        }
    } else if (node.is_object()) {
        for (auto it = node.begin(); it != node.end(); ++it) {
            const std::string& key = it.key();
            const json& el = it.value();
            std::string name = key;
            if (el.is_object() && el.contains("name") && el["name"].is_string()) {
                name = el["name"].get<std::string>();
            }
            cb(name, el);
        }
    } else {
        throw std::runtime_error(
            "LibraryContext: '" + kind + "' must be a JSON object or array");
    }
}

} // namespace

void LibraryContext::loadFromJson(const json& data, LoadMode mode) {
    _mode = mode;

    if (data.contains("coreShapes")) {
        std::map<std::string, MAS::CoreShape> out;
        forEachEntry(data["coreShapes"], "coreShapes", [&](const std::string& name, const json& jf) {
            MAS::CoreShape shape;
            json normalized = jf;
            if (normalized.contains("family") && normalized["family"].is_string()) {
                std::string family = normalized["family"];
                std::transform(family.begin(), family.end(), family.begin(), ::toupper);
                std::replace(family.begin(), family.end(), ' ', '_');
                normalized["family"] = family;
            }
            from_json(normalized, shape);
            out[name] = shape;
        });
        _coreShapes = std::move(out);
    }

    if (data.contains("coreMaterials")) {
        std::map<std::string, MAS::CoreMaterial> out;
        forEachEntry(data["coreMaterials"], "coreMaterials", [&](const std::string& name, const json& jf) {
            try {
                MAS::CoreMaterial m(jf);
                out[name] = m;
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "LibraryContext: failed to parse coreMaterial '" + name + "': " + e.what());
            }
        });
        _coreMaterials = std::move(out);
    }

    if (data.contains("cores")) {
        std::vector<Core> out;
        forEachEntry(data["cores"], "cores", [&](const std::string& name, const json& jf) {
            try {
                out.emplace_back(jf);
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "LibraryContext: failed to parse core '" + name + "': " + e.what());
            }
        });
        _cores = std::move(out);
    }

    if (data.contains("wires")) {
        std::map<std::string, Wire> out;
        forEachEntry(data["wires"], "wires", [&](const std::string& name, const json& jf) {
            try {
                Wire w(jf);
                out[name] = w;
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "LibraryContext: failed to parse wire '" + name + "': " + e.what());
            }
        });
        _wires = std::move(out);
    }

    if (data.contains("bobbins")) {
        std::map<std::string, Bobbin> out;
        forEachEntry(data["bobbins"], "bobbins", [&](const std::string& name, const json& jf) {
            try {
                Bobbin b(jf);
                out[name] = b;
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "LibraryContext: failed to parse bobbin '" + name + "': " + e.what());
            }
        });
        _bobbins = std::move(out);
    }

    if (data.contains("insulationMaterials")) {
        std::map<std::string, InsulationMaterial> out;
        forEachEntry(data["insulationMaterials"], "insulationMaterials", [&](const std::string& name, const json& jf) {
            try {
                InsulationMaterial im(jf);
                out[name] = im;
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "LibraryContext: failed to parse insulationMaterial '" + name + "': " + e.what());
            }
        });
        _insulationMaterials = std::move(out);
    }

    if (data.contains("wireMaterials")) {
        std::map<std::string, MAS::WireMaterial> out;
        forEachEntry(data["wireMaterials"], "wireMaterials", [&](const std::string& name, const json& jf) {
            try {
                MAS::WireMaterial wm(jf);
                out[name] = wm;
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "LibraryContext: failed to parse wireMaterial '" + name + "': " + e.what());
            }
        });
        _wireMaterials = std::move(out);
    }
}

void LibraryContext::loadFromString(const std::string& jsonText, LoadMode mode) {
    json data;
    try {
        data = json::parse(jsonText);
    } catch (const std::exception& e) {
        throw std::runtime_error(
            std::string("LibraryContext: invalid JSON in loadFromString: ") + e.what());
    }
    loadFromJson(data, mode);
}

void LibraryContext::loadFromFile(const std::string& path, LoadMode mode) {
    std::ifstream f(path);
    if (!f.is_open()) {
        throw std::runtime_error("LibraryContext: cannot open file: " + path);
    }
    std::ostringstream ss;
    ss << f.rdbuf();
    loadFromString(ss.str(), mode);
}

void LibraryContext::clear() {
    _cores.reset();
    _coreMaterials.reset();
    _coreShapes.reset();
    _wires.reset();
    _bobbins.reset();
    _insulationMaterials.reset();
    _wireMaterials.reset();
}

bool LibraryContext::empty() const {
    return !_cores && !_coreMaterials && !_coreShapes && !_wires
        && !_bobbins && !_insulationMaterials && !_wireMaterials;
}

// --- Scope (RAII swap) -----------------------------------------------------

namespace {

// Recompute coreShapeFamiliesInDatabase from the current coreShapeDatabase.
void rebuildShapeFamilies() {
    coreShapeFamiliesInDatabase.clear();
    for (const auto& [name, shape] : coreShapeDatabase) {
        auto fam = shape.get_family();
        if (std::find(coreShapeFamiliesInDatabase.begin(),
                      coreShapeFamiliesInDatabase.end(),
                      fam) == coreShapeFamiliesInDatabase.end()) {
            coreShapeFamiliesInDatabase.push_back(fam);
        }
    }
}

} // namespace

LibraryContext::Scope::Scope(const LibraryContext* ctx) {
    if (!ctx || ctx->empty()) {
        return;  // nothing to do; global state untouched
    }
    _active = true;
    _coreBackup = coreDatabase;
    _coreMaterialBackup = coreMaterialDatabase;
    _coreShapeBackup = coreShapeDatabase;
    _coreShapeFamilyBackup = coreShapeFamiliesInDatabase;
    _wireBackup = wireDatabase;
    _bobbinBackup = bobbinDatabase;
    _insulationMaterialBackup = insulationMaterialDatabase;
    _wireMaterialBackup = wireMaterialDatabase;

    const bool replace = (ctx->mode() == LoadMode::Replace);

    if (ctx->cores()) {
        if (replace) coreDatabase.clear();
        for (const auto& c : *ctx->cores()) coreDatabase.push_back(c);
    }
    if (ctx->coreMaterials()) {
        if (replace) coreMaterialDatabase.clear();
        for (const auto& [k, v] : *ctx->coreMaterials()) coreMaterialDatabase[k] = v;
    }
    if (ctx->coreShapes()) {
        if (replace) coreShapeDatabase.clear();
        for (const auto& [k, v] : *ctx->coreShapes()) coreShapeDatabase[k] = v;
        rebuildShapeFamilies();
    }
    if (ctx->wires()) {
        if (replace) wireDatabase.clear();
        for (const auto& [k, v] : *ctx->wires()) wireDatabase[k] = v;
    }
    if (ctx->bobbins()) {
        if (replace) bobbinDatabase.clear();
        for (const auto& [k, v] : *ctx->bobbins()) bobbinDatabase[k] = v;
    }
    if (ctx->insulationMaterials()) {
        if (replace) insulationMaterialDatabase.clear();
        for (const auto& [k, v] : *ctx->insulationMaterials()) insulationMaterialDatabase[k] = v;
    }
    if (ctx->wireMaterials()) {
        if (replace) wireMaterialDatabase.clear();
        for (const auto& [k, v] : *ctx->wireMaterials()) wireMaterialDatabase[k] = v;
    }
}

LibraryContext::Scope::~Scope() {
    if (!_active) return;
    coreDatabase = std::move(_coreBackup);
    coreMaterialDatabase = std::move(_coreMaterialBackup);
    coreShapeDatabase = std::move(_coreShapeBackup);
    coreShapeFamiliesInDatabase = std::move(_coreShapeFamilyBackup);
    wireDatabase = std::move(_wireBackup);
    bobbinDatabase = std::move(_bobbinBackup);
    insulationMaterialDatabase = std::move(_insulationMaterialBackup);
    wireMaterialDatabase = std::move(_wireMaterialBackup);
}

LibraryContext::Scope::Scope(Scope&& other) noexcept
    : _active(other._active),
      _coreBackup(std::move(other._coreBackup)),
      _coreMaterialBackup(std::move(other._coreMaterialBackup)),
      _coreShapeBackup(std::move(other._coreShapeBackup)),
      _coreShapeFamilyBackup(std::move(other._coreShapeFamilyBackup)),
      _wireBackup(std::move(other._wireBackup)),
      _bobbinBackup(std::move(other._bobbinBackup)),
      _insulationMaterialBackup(std::move(other._insulationMaterialBackup)),
      _wireMaterialBackup(std::move(other._wireMaterialBackup)) {
    other._active = false;
}

LibraryContext::Scope& LibraryContext::Scope::operator=(Scope&& other) noexcept {
    if (this != &other) {
        this->~Scope();
        new (this) Scope(std::move(other));
    }
    return *this;
}

// --- Type-filter predicates ------------------------------------------------

namespace {

std::string toUpperCopy(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    return s;
}

template <typename E>
bool acceptsEnum(const TypeFilterSet& f, E value) {
    if (f.empty()) return true;
    std::string name = std::string(magic_enum::enum_name(value));
    std::string upper = toUpperCopy(name);

    auto contains = [&](const std::set<std::string>& s) {
        for (const auto& v : s) {
            if (toUpperCopy(v) == upper) return true;
        }
        return false;
    };

    if (!f.allowed.empty() && !contains(f.allowed)) return false;
    if (!f.blocked.empty() && contains(f.blocked)) return false;
    return true;
}

} // namespace

bool acceptsCoreShapeFamily(const TypeFilterSet& f, MAS::CoreShapeFamily family) {
    return acceptsEnum(f, family);
}

bool acceptsCoreMaterialType(const TypeFilterSet& f, MAS::MaterialType type) {
    return acceptsEnum(f, type);
}

bool acceptsWireType(const TypeFilterSet& f, MAS::WireType type) {
    return acceptsEnum(f, type);
}

std::vector<Core> filterCoresByConstraints(const std::vector<Core>& cores,
                                            const AdviserConstraints& constraints) {
    if (constraints.shapeFamily.empty() && constraints.coreMaterialType.empty()) {
        return cores;
    }
    std::vector<Core> out;
    out.reserve(cores.size());
    for (const auto& core : cores) {
        const auto& shape = core.get_functional_description().get_shape();
        bool shapeOk = true;
        if (!constraints.shapeFamily.empty()) {
            if (std::holds_alternative<MAS::CoreShape>(shape)) {
                shapeOk = acceptsCoreShapeFamily(constraints.shapeFamily,
                                                  std::get<MAS::CoreShape>(shape).get_family());
            } else {
                // Shape is referenced by name only; look up in the active database.
                std::string name = std::get<std::string>(shape);
                auto it = coreShapeDatabase.find(name);
                if (it == coreShapeDatabase.end()) {
                    shapeOk = false;
                } else {
                    shapeOk = acceptsCoreShapeFamily(constraints.shapeFamily, it->second.get_family());
                }
            }
        }
        if (!shapeOk) continue;

        bool materialOk = true;
        if (!constraints.coreMaterialType.empty()) {
            const auto& mat = core.get_functional_description().get_material();
            MAS::CoreMaterial resolved;
            bool found = false;
            if (std::holds_alternative<MAS::CoreMaterial>(mat)) {
                resolved = std::get<MAS::CoreMaterial>(mat);
                found = true;
            } else {
                std::string name = std::get<std::string>(mat);
                auto it = coreMaterialDatabase.find(name);
                if (it != coreMaterialDatabase.end()) {
                    resolved = it->second;
                    found = true;
                }
            }
            if (!found) {
                materialOk = false;
            } else {
                materialOk = acceptsCoreMaterialType(constraints.coreMaterialType,
                                                     resolved.get_material());
            }
        }
        if (!materialOk) continue;

        out.push_back(core);
    }
    return out;
}

std::map<std::string, Wire> filterWiresByConstraints(
    const std::map<std::string, Wire>& wires,
    const AdviserConstraints& constraints) {
    if (constraints.wireType.empty()) return wires;
    std::map<std::string, Wire> out;
    for (const auto& [name, wire] : wires) {
        if (acceptsWireType(constraints.wireType, wire.get_type())) {
            out.emplace(name, wire);
        }
    }
    return out;
}

} // namespace OpenMagnetics
