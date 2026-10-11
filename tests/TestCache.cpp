#include "constructive_models/Mas.h"
#include "support/Utils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using json = nlohmann::json;
#include <typeinfo>

using namespace MAS;
using namespace OpenMagnetics;

namespace {
TEST_CASE("Magnetic_Cache_Load", "[support][cache][smoke-test]") {
    std::string masString = R"({"outputs": [], "inputs": {"designRequirements": {"isolationSides": ["primary" ], "magnetizingInductance": {"nominal": 0.00039999999999999996 }, "name": "My Design Requirements", "turnsRatios": [{"nominal": 1} ] }, "operatingPoints": [{"conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"frequency": 100000, "current": {"processed": {"label": "triangular", "peakToPeak": 0.5, "offset": 0, "dutyCycle": 0.5 } }, "voltage": {"processed": {"label": "rectangular", "peakToPeak": 20, "offset": 0, "dutyCycle": 0.5 } } } ], "name": "Operating Point No. 1" } ] }, "magnetic": {"coil": {"bobbin": "basic", "functionalDescription":[{"name": "Primary", "numberTurns": 4, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 4, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] }, "core": {"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "twoPieceSet", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } }, "manufacturerInfo": {"name": "", "reference": "Example" } } })";
    json masJson = json::parse(masString);
    OpenMagnetics::Mas mas(masJson);

    auto magnetic = mas.get_magnetic();

    magneticsCache.load("A", magnetic);
    REQUIRE(magneticsCache.size() == 1);
    REQUIRE(magneticsCache.energy_cache_size() == 0);
}

TEST_CASE("Magnetic_Cache_Compute_Energy", "[support][cache][smoke-test]") {
    {
        std::string masString = R"({"outputs": [], "inputs": {"designRequirements": {"isolationSides": ["primary" ], "magnetizingInductance": {"nominal": 0.00039999999999999996 }, "name": "My Design Requirements", "turnsRatios": [{"nominal": 1} ] }, "operatingPoints": [{"conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"frequency": 100000, "current": {"processed": {"label": "triangular", "peakToPeak": 0.5, "offset": 0, "dutyCycle": 0.5 } }, "voltage": {"processed": {"label": "rectangular", "peakToPeak": 20, "offset": 0, "dutyCycle": 0.5 } } } ], "name": "Operating Point No. 1" } ] }, "magnetic": {"coil": {"bobbin": "basic", "functionalDescription":[{"name": "Primary", "numberTurns": 4, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 4, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] }, "core": {"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "twoPieceSet", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } }, "manufacturerInfo": {"name": "", "reference": "Example" } } })";
        json masJson = json::parse(masString);
        OpenMagnetics::Mas mas(masJson);
        auto magnetic = mas.get_magnetic();
        magneticsCache.load("A", magnetic);
    }
    {
        std::string masString = R"({"outputs": [], "inputs": {"designRequirements": {"isolationSides": ["primary" ], "magnetizingInductance": {"nominal": 0.00039999999999999996 }, "name": "My Design Requirements", "turnsRatios": [{"nominal": 1} ] }, "operatingPoints": [{"conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"frequency": 100000, "current": {"processed": {"label": "triangular", "peakToPeak": 0.5, "offset": 0, "dutyCycle": 0.5 } }, "voltage": {"processed": {"label": "rectangular", "peakToPeak": 20, "offset": 0, "dutyCycle": 0.5 } } } ], "name": "Operating Point No. 1" } ] }, "magnetic": {"coil": {"bobbin": "basic", "functionalDescription":[{"name": "Primary", "numberTurns": 4, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 4, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] }, "core": {"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "twoPieceSet", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.00001 }], "numberStacks": 1 } }, "manufacturerInfo": {"name": "", "reference": "Example" } } })";
        json masJson = json::parse(masString);
        OpenMagnetics::Mas mas(masJson);
        auto magnetic = mas.get_magnetic();
        magneticsCache.load("A with different gap", magnetic);
    }
    REQUIRE(magneticsCache.size() == 2);

    magneticsCache.autocomplete_magnetics();

    auto references = magneticsCache.references();
    REQUIRE(references[0] == "A");
    REQUIRE(references[1] == "A with different gap");
    REQUIRE(magneticsCache.energy_cache_size() == 0);

    magneticsCache.compute_energy_cache();
    REQUIRE(magneticsCache.energy_cache_size() == 2);
}

// ABT #1496: magnetic_autocomplete throws on a winding that does not fit its window. A catalogue
// load must neither keep such a record wound past its window nor drop it quietly, and stopping at
// the first one would hide every later one: the load fails ONCE, naming every unfit record, and
// leaves the cache as it was.
TEST_CASE("Magnetic_Cache_Autocomplete_Names_Every_Unfit_Record", "[support][cache][abt1496]") {
    magneticsCache.clear();
    auto magneticWithTurns = [](int64_t numberTurns) {
        std::string masString = R"({"outputs": [], "inputs": {"designRequirements": {"isolationSides": ["primary" ], "magnetizingInductance": {"nominal": 0.0004 }, "name": "My Design Requirements", "turnsRatios": [{"nominal": 1} ] }, "operatingPoints": [{"conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"frequency": 100000, "current": {"processed": {"label": "triangular", "peakToPeak": 0.5, "offset": 0, "dutyCycle": 0.5 } }, "voltage": {"processed": {"label": "rectangular", "peakToPeak": 20, "offset": 0, "dutyCycle": 0.5 } } } ], "name": "Operating Point No. 1" } ] }, "magnetic": {"coil": {"bobbin": "basic", "functionalDescription":[{"name": "Primary", "numberTurns": 4, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 4, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] }, "core": {"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "twoPieceSet", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } }, "manufacturerInfo": {"name": "", "reference": "Example" } } })";
        json masJson = json::parse(masString);
        masJson["magnetic"]["coil"]["functionalDescription"][0]["numberTurns"] = numberTurns;
        masJson["magnetic"]["coil"]["functionalDescription"][1]["numberTurns"] = numberTurns;
        OpenMagnetics::Mas mas(masJson);
        return mas.get_magnetic();
    };
    // 4 + 4 turns of a 1 mm wire fit a PQ 32/20; 400 + 400 do not.
    magneticsCache.load("fits", magneticWithTurns(4));
    magneticsCache.load("unfit one", magneticWithTurns(400));
    magneticsCache.load("unfit two", magneticWithTurns(400));
    REQUIRE(magneticsCache.size() == 3);

    bool threw = false;
    try {
        magneticsCache.autocomplete_magnetics();
    }
    catch (const OpenMagnetics::OpenMagneticsException& e) {
        threw = true;
        const std::string what = e.what();
        INFO(what);
        CHECK(e.code() == OpenMagnetics::ErrorCode::COIL_WINDING_ERROR);
        CHECK(what.find("2 of 3 records") != std::string::npos);
        CHECK(what.find("'unfit one'") != std::string::npos);
        CHECK(what.find("'unfit two'") != std::string::npos);
        CHECK(what.find("'fits'") == std::string::npos);
    }
    CHECK(threw);
    // Nothing was replaced: every record is still the one that was loaded, unwound.
    REQUIRE(magneticsCache.size() == 3);
    for (const auto& reference : magneticsCache.references()) {
        CHECK_FALSE(magneticsCache.read(reference).get_coil().get_turns_description());
    }
    magneticsCache.clear();
}
}  // namespace
