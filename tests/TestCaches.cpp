#include "constructive_models/Mas.h"
#include "support/Caches.h"

#include <UnitTest++.h>

using json = nlohmann::json;
#include <typeinfo>

using namespace MAS;
using namespace OpenMagnetics;

SUITE(Caches) {
    TEST(Magnetic_Cache_Load) {
        std::string masString = R"({"outputs": [], "inputs": {"designRequirements": {"isolationSides": ["primary" ], "magnetizingInductance": {"nominal": 0.00039999999999999996 }, "name": "My Design Requirements", "turnsRatios": [{"nominal": 1} ] }, "operatingPoints": [{"conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"frequency": 100000, "current": {"processed": {"label": "Triangular", "peakToPeak": 0.5, "offset": 0, "dutyCycle": 0.5 } }, "voltage": {"processed": {"label": "Rectangular", "peakToPeak": 20, "offset": 0, "dutyCycle": 0.5 } } } ], "name": "Operating Point No. 1" } ] }, "magnetic": {"coil": {"bobbin": "Basic", "functionalDescription":[{"name": "Primary", "numberTurns": 4, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 4, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] }, "core": {"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } }, "manufacturerInfo": {"name": "", "reference": "Example" } } })";
        json masJson = json::parse(masString);
        OpenMagnetics::Mas mas(masJson);

        auto magnetic = mas.get_magnetic();

        load_magnetic_in_cache("A", magnetic);
        CHECK(magneticsCache.size() == 1);
        CHECK(magneticEnergyCache.size() == 0);
    }

    TEST(Magnetic_Cache_Compute_Energy) {
        {
            std::string masString = R"({"outputs": [], "inputs": {"designRequirements": {"isolationSides": ["primary" ], "magnetizingInductance": {"nominal": 0.00039999999999999996 }, "name": "My Design Requirements", "turnsRatios": [{"nominal": 1} ] }, "operatingPoints": [{"conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"frequency": 100000, "current": {"processed": {"label": "Triangular", "peakToPeak": 0.5, "offset": 0, "dutyCycle": 0.5 } }, "voltage": {"processed": {"label": "Rectangular", "peakToPeak": 20, "offset": 0, "dutyCycle": 0.5 } } } ], "name": "Operating Point No. 1" } ] }, "magnetic": {"coil": {"bobbin": "Basic", "functionalDescription":[{"name": "Primary", "numberTurns": 4, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 4, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] }, "core": {"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } }, "manufacturerInfo": {"name": "", "reference": "Example" } } })";
            json masJson = json::parse(masString);
            OpenMagnetics::Mas mas(masJson);
            auto magnetic = mas.get_magnetic();
            load_magnetic_in_cache("A", magnetic);
        }
        {
            std::string masString = R"({"outputs": [], "inputs": {"designRequirements": {"isolationSides": ["primary" ], "magnetizingInductance": {"nominal": 0.00039999999999999996 }, "name": "My Design Requirements", "turnsRatios": [{"nominal": 1} ] }, "operatingPoints": [{"conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"frequency": 100000, "current": {"processed": {"label": "Triangular", "peakToPeak": 0.5, "offset": 0, "dutyCycle": 0.5 } }, "voltage": {"processed": {"label": "Rectangular", "peakToPeak": 20, "offset": 0, "dutyCycle": 0.5 } } } ], "name": "Operating Point No. 1" } ] }, "magnetic": {"coil": {"bobbin": "Basic", "functionalDescription":[{"name": "Primary", "numberTurns": 4, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 4, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] }, "core": {"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.00001 }], "numberStacks": 1 } }, "manufacturerInfo": {"name": "", "reference": "Example" } } })";
            json masJson = json::parse(masString);
            OpenMagnetics::Mas mas(masJson);
            auto magnetic = mas.get_magnetic();
            load_magnetic_in_cache("A with different gap", magnetic);
        }
        CHECK(magneticsCache.size() == 2);

        autocomplete_magnetics_in_cache();

        auto references = get_magnetic_cache_references();
        CHECK(references[0] == "A");
        CHECK(references[1] == "A with different gap");
        CHECK(magneticEnergyCache.size() == 0);

        compute_energy_cache();
        CHECK(magneticEnergyCache.size() == 2);
    }
}