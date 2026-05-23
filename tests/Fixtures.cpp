#include "Fixtures.h"

#include "constructive_models/Bobbin.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/MasMigration.h"
#include "physical_models/MagnetizingInductance.h"
#include "processors/Inputs.h"
#include "support/Utils.h"

#include <algorithm>
#include <atomic>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <sstream>
#include <stdexcept>

namespace OpenMagneticsTesting::fixtures {

namespace {

using nlohmann::json;
using std::filesystem::path;

// In-memory registry. Built once, then never mutated except via
// reload_for_testing() (which holds the same mutex). We use
// atomic<bool>+mutex rather than std::once_flag because reload_for_testing()
// must be able to reset the "already initialised" signal, and once_flag is
// neither copyable nor assignable.
std::map<std::string, json> g_registry;
std::atomic<bool> g_initialized{false};
std::mutex g_initMutex;

// Locate the fixtures directory. The convention is
// tests/fixtures/ relative to the source tree. We walk a few
// candidate roots so the loader works from both CTest (build dir)
// and direct invocation (source dir).
path locate_fixtures_dir() {
    // Honor an env var override first — useful for sandboxed CI.
    if (const char* env = std::getenv("MKF_FIXTURES_DIR")) {
        path p{env};
        if (std::filesystem::is_directory(p)) return p;
        throw std::runtime_error(
            "MKF_FIXTURES_DIR is set but doesn't point at a directory: " +
            p.string());
    }

    // Probe candidates that match where tests typically run.
    const std::vector<path> candidates = {
        path{"tests"} / "fixtures",
        path{".."} / "tests" / "fixtures",
        path{"../.."} / "tests" / "fixtures",
    };
    for (const auto& c : candidates) {
        if (std::filesystem::is_directory(c)) return c;
    }

    // Last-ditch: walk upward from CWD looking for a `tests/fixtures`.
    path cwd = std::filesystem::current_path();
    for (int i = 0; i < 5; ++i) {
        path candidate = cwd / "tests" / "fixtures";
        if (std::filesystem::is_directory(candidate)) return candidate;
        if (!cwd.has_parent_path()) break;
        cwd = cwd.parent_path();
    }

    throw std::runtime_error(
        "Could not locate tests/fixtures/ directory. Set the "
        "MKF_FIXTURES_DIR env var to its absolute path.");
}

// Repair a known ill-formed UTF-8 pattern (lone 0xB5 Latin-1 µ ->
// proper 0xC2 0xB5). Mirrors the in-mas_loader fix so on-disk
// fixtures can be edited with non-UTF-8 tooling and still load.
void repair_utf8(std::string& s) {
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (static_cast<unsigned char>(s[i]) == 0xB5 &&
            (i == 0 || static_cast<unsigned char>(s[i - 1]) != 0xC2)) {
            s.insert(i, 1, static_cast<char>(0xC2));
            ++i;
        }
    }
}

void load_one_file(const path& filePath) {
    std::ifstream f(filePath);
    if (!f) {
        throw std::runtime_error(
            "Could not open fixture file: " + filePath.string());
    }

    std::string line;
    std::size_t lineNo = 0;
    while (std::getline(f, line)) {
        ++lineNo;
        // Skip blank and comment-y lines so users can annotate ndjson
        // sources lightly (note: real JSON doesn't have comments, so
        // any '#' or '//' prefix line is silently dropped).
        if (line.empty()) continue;
        // Trim leading whitespace and check for comment markers.
        const auto firstNonWs = line.find_first_not_of(" \t");
        if (firstNonWs == std::string::npos) continue;
        if (line[firstNonWs] == '#' ||
            (line[firstNonWs] == '/' && firstNonWs + 1 < line.size() &&
             line[firstNonWs + 1] == '/')) {
            continue;
        }

        repair_utf8(line);

        json parsed;
        try {
            parsed = json::parse(line);
        } catch (const std::exception& e) {
            throw std::runtime_error(
                "Fixture parse error at " + filePath.string() + ":" +
                std::to_string(lineNo) + ": " + e.what());
        }

        if (!parsed.is_object()) {
            throw std::runtime_error(
                "Fixture must be a JSON object at " + filePath.string() +
                ":" + std::to_string(lineNo) +
                " (got " + std::string(parsed.type_name()) + ")");
        }
        if (!parsed.contains("name") || !parsed["name"].is_string()) {
            throw std::runtime_error(
                "Fixture missing string \"name\" at " + filePath.string() +
                ":" + std::to_string(lineNo));
        }
        if (!parsed.contains("data")) {
            throw std::runtime_error(
                "Fixture missing \"data\" at " + filePath.string() + ":" +
                std::to_string(lineNo));
        }

        std::string name = parsed["name"].get<std::string>();
        if (g_registry.count(name) > 0) {
            throw std::runtime_error(
                "Duplicate fixture name \"" + name +
                "\" — second occurrence at " + filePath.string() + ":" +
                std::to_string(lineNo));
        }
        g_registry.emplace(std::move(name), std::move(parsed["data"]));
    }
}

void init_registry() {
    const path root = locate_fixtures_dir();

    // Collect candidate files first so iteration order is stable —
    // makes duplicate-name errors deterministic.
    std::vector<path> files;
    for (const auto& entry :
            std::filesystem::recursive_directory_iterator(root)) {
        if (!entry.is_regular_file()) continue;
        if (entry.path().extension() == ".ndjson") files.push_back(entry.path());
    }
    std::sort(files.begin(), files.end());

    for (const auto& f : files) {
        load_one_file(f);
    }
}

void ensure_initialized() {
    // Double-checked locking on the atomic flag — cheap fast path, only
    // takes the mutex when we genuinely need to populate.
    if (g_initialized.load(std::memory_order_acquire)) return;
    std::lock_guard lk(g_initMutex);
    if (g_initialized.load(std::memory_order_relaxed)) return;
    init_registry();
    g_initialized.store(true, std::memory_order_release);
}

// Same three-tier inductance recovery as TestingUtils::mas_loader.
// Kept duplicated for now (rather than refactoring mas_loader to
// share) so this code can land as a self-contained addition without
// touching the existing loader's call-sites. A follow-up commit can
// fold the two implementations together.
OpenMagnetics::Mas mas_from_payload(const json& payload, const std::string& fixtureName) {
    auto masJson = payload;
    OpenMagnetics::compat::migrate_pre_1_0(masJson);

    auto inputsJson = masJson["inputs"];
    auto magneticJson = masJson["magnetic"];
    auto outputsJson = masJson["outputs"];

    auto magnetic = OpenMagnetics::Magnetic(magneticJson);
    if (magneticJson["coil"]["bobbin"] == "basic") {
        auto bobbinData = OpenMagnetics::Bobbin::create_quick_bobbin(
            magnetic.get_mutable_core(), false);
        to_json(magneticJson["coil"]["bobbin"], bobbinData);
    }
    auto coil = OpenMagnetics::Coil(magneticJson["coil"]);

    std::vector<OpenMagnetics::Outputs> outputs;
    if (outputsJson != nullptr) {
        outputs = outputsJson.get<std::vector<OpenMagnetics::Outputs>>();
    }

    OpenMagnetics::Inputs inputs;
    std::vector<double> magnetizingInductancePerPoint;
    for (const auto& output : outputs) {
        if (output.get_inductance()) {
            auto mi = OpenMagnetics::resolve_dimensional_values(
                output.get_inductance()->get_magnetizing_inductance()
                    .get_magnetizing_inductance());
            magnetizingInductancePerPoint.push_back(mi);
        }
    }

    if (!magnetizingInductancePerPoint.empty()) {
        inputs = OpenMagnetics::Inputs(inputsJson, true,
                                        magnetizingInductancePerPoint);
    } else {
        try {
            OpenMagnetics::MagnetizingInductance model;
            double L = model.calculate_inductance_from_number_turns_and_gapping(
                            magnetic.get_core(), magnetic.get_coil())
                            .get_magnetizing_inductance()
                            .get_nominal()
                            .value();
            inputs = OpenMagnetics::Inputs(inputsJson, true, L);
        } catch (const std::exception& e) {
            std::cerr << "[fixtures] " << fixtureName
                      << ": couldn't recover magnetizing inductance ("
                      << e.what()
                      << ") — falling back to default-construct Inputs.\n";
            inputs = OpenMagnetics::Inputs(inputsJson, true);
        }
    }

    OpenMagnetics::Mas mas;
    mas.set_inputs(inputs);
    mas.set_magnetic(magnetic);
    mas.set_outputs(outputs);
    return mas;
}

// Compose a "did you mean…" suggestion list when a fixture isn't
// found. Cheap O(N) scan over the registry — only triggers on the
// error path so this never impacts steady-state lookups.
std::string suggest_close_matches(const std::string& name) {
    std::vector<std::string> suggestions;
    for (const auto& [key, _] : g_registry) {
        // Trivial substring heuristic — keeps the loader dependency
        // graph tiny (no rapidfuzz here). Good enough for typos.
        if (key.find(name) != std::string::npos ||
            name.find(key) != std::string::npos) {
            suggestions.push_back(key);
            if (suggestions.size() >= 5) break;
        }
    }
    if (suggestions.empty()) return "";
    std::ostringstream os;
    os << " (did you mean:";
    for (const auto& s : suggestions) os << " \"" << s << "\"";
    os << "?)";
    return os.str();
}

}  // namespace

json get_json(const std::string& name) {
    ensure_initialized();
    auto it = g_registry.find(name);
    if (it == g_registry.end()) {
        throw std::runtime_error(
            "Unknown fixture: \"" + name + "\"" + suggest_close_matches(name) +
            " — registry has " + std::to_string(g_registry.size()) +
            " entries");
    }
    return it->second;  // copy on purpose — fixtures are immutable, callers may mutate.
}

OpenMagnetics::Mas get_mas(const std::string& name) {
    return mas_from_payload(get_json(name), name);
}

std::vector<std::string> list_names() {
    ensure_initialized();
    std::vector<std::string> names;
    names.reserve(g_registry.size());
    for (const auto& [k, _] : g_registry) names.push_back(k);
    return names;
}

void reload_for_testing() {
    std::lock_guard lk(g_initMutex);
    g_registry.clear();
    g_initialized.store(false, std::memory_order_release);
    init_registry();
    g_initialized.store(true, std::memory_order_release);
}

}  // namespace OpenMagneticsTesting::fixtures
