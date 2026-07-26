// Validates sourced shape geometries for the new families (EER, EP-derivatives) by building
// each shape record and comparing MKF-computed effective Ae/Le against the datasheet values.
// Data files live in src/tools/vendor_cores/ (produced by the sourcing agents):
//   <fam>_shapes.ndjson  — MAS core-shape records (dims in metres)
//   <fam>_validation.csv  — name,datasheet_ae_mm2,datasheet_le_mm[,...]
// A correct A-F -> MKF-convention mapping yields Ae within ~15% and Le within ~12% of the
// datasheet (mechanical centre-leg vs effective area differ a few %). A gross mismatch (>~25%)
// means the dimensions were mislabelled.
#include "constructive_models/Core.h"
#include "json.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <fstream>
#include <sstream>
#include <map>

using namespace MAS;
using namespace OpenMagnetics;
using json = nlohmann::json;

static const std::string DIR = std::string(getenv("VENDOR_CORES_DIR") ?
    getenv("VENDOR_CORES_DIR") : "src/tools/vendor_cores");

static std::map<std::string, std::pair<double,double>> load_validation(const std::string& path) {
    std::map<std::string, std::pair<double,double>> out;
    std::ifstream f(path);
    std::string line; bool header = true;
    while (std::getline(f, line)) {
        if (header) { header = false; continue; }
        if (line.empty()) continue;
        std::stringstream ss(line); std::string name, ae, le;
        std::getline(ss, name, ','); std::getline(ss, ae, ','); std::getline(ss, le, ',');
        try { out[name] = {std::stod(ae), std::stod(le)}; } catch (...) {}
    }
    return out;
}

static void validate_family(const std::string& shapesFile, const std::string& valFile) {
    std::ifstream sf(DIR + "/" + shapesFile);
    if (!sf.good()) { WARN("skip: " << shapesFile << " not present"); return; }
    auto expected = load_validation(DIR + "/" + valFile);
    std::string line;
    int checked = 0;
    while (std::getline(sf, line)) {
        if (line.empty()) continue;
        json js = json::parse(line);
        std::string name = js.at("name");
        // Build the FULL assembled core (two-piece set) so effective params match the
        // datasheet's assembled-core Ae/Le — the per-piece partial params would report ~half Le.
        json coreJson = {
            {"name", name},
            {"functionalDescription", {
                {"type", "two-piece set"},
                {"material", "3C90"},
                {"shape", js},
                {"gapping", json::array()},
                {"numberStacks", 1}
            }}
        };
        Core core(coreJson, true);
        auto ep = core.get_processed_description()->get_effective_parameters();
        double aeMm2 = ep.get_effective_area() * 1e6;
        double leMm  = ep.get_effective_length() * 1e3;
        INFO("shape=" << name << " MKF Ae=" << aeMm2 << " Le=" << leMm);
        REQUIRE(aeMm2 > 0.0);
        REQUIRE(leMm > 0.0);
        auto it = expected.find(name);
        if (it != expected.end()) {
            double eAe = it->second.first, eLe = it->second.second;
            INFO("datasheet Ae=" << eAe << " Le=" << eLe);
            REQUIRE_THAT(aeMm2, Catch::Matchers::WithinRel(eAe, 0.20));
            REQUIRE_THAT(leMm,  Catch::Matchers::WithinRel(eLe, 0.15));
            checked++;
        }
    }
    INFO("validated " << checked << " shapes in " << shapesFile);
    CHECK(checked > 0);
}

TEST_CASE("EER sourced geometry matches datasheet Ae/Le", "[newgeom][eer]") {
    validate_family("eer_shapes.ndjson", "eer_validation.csv");
}
TEST_CASE("EP-derivative sourced geometry matches datasheet Ae/Le", "[newgeom][epderiv]") {
    validate_family("epderiv_shapes.ndjson", "epderiv_validation.csv");
}
