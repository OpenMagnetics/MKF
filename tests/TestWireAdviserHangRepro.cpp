/**
 * Repro for the Wire Configuration → Advise All hang.
 *
 * Fixture: tests/payloads/wire_adviser_hang.json
 *   Captured from the frontend Magnetic Builder before it calls
 *   mkf.calculate_advised_coil(JSON.stringify(mas)).
 *
 * Setup mirrors WebLibMKF/src/libMKF.cpp:2860 calculate_advised_coil:
 *   - Force every winding's wire to "Dummy" and clear bobbin/coil descriptions
 *     (the frontend strips those when calling Wire Adviser).
 *   - Run CoilAdviser::get_advised_coil(mas, 1).
 *   - Watchdog at 60 s so a true hang doesn't block the test suite.
 *
 * Goal: surface the slow loop (file:line / iteration count) so the user can
 * decide on a fix. The REQUIRE is the same contract as the frontend button:
 * the adviser MUST return at least one design within the watchdog window.
 */
#include "advisers/CoilAdviser.h"
#include "processors/Inputs.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <source_location>
#include <sstream>

using namespace MAS;
using namespace OpenMagnetics;

TEST_CASE("Test_CoilAdviser_WireAdviser_Hang_Repro",
          "[adviser][coil-adviser][wire-adviser][hang][bug-repro]") {
    clear_databases();

    // Load the captured payload (MAS-shaped JSON).
    auto path = std::filesystem::path{std::source_location::current().file_name()}
                    .parent_path()
                    .append("payloads")
                    .append("wire_adviser_hang.json");
    std::ifstream f(path);
    REQUIRE(f.is_open());
    std::stringstream buf;
    buf << f.rdbuf();
    auto raw = buf.str();

    json masJson = json::parse(raw);

    // Mirror libMKF::calculate_advised_coil(...): wires set to "Dummy",
    // clear coil descriptions, clear bobbin/coil derived fields.
    auto& functional = masJson["magnetic"]["coil"]["functionalDescription"];
    for (auto& winding : functional) {
        winding["wire"] = "Dummy";
    }
    masJson["magnetic"]["coil"]["sectionsDescription"] = nullptr;
    masJson["magnetic"]["coil"]["layersDescription"]   = nullptr;
    masJson["magnetic"]["coil"]["turnsDescription"]    = nullptr;

    OpenMagnetics::Mas mas(masJson);

    // Settings that match the WebLibMKF Wire Adviser entry path. The frontend
    // doesn't explicitly tweak these for calculate_advised_coil — defaults
    // apply. Force margin tape ON (the broader space the adviser explores).
    auto& settings = Settings::GetInstance();

    CoilAdviser coilAdviser;

    auto start = std::chrono::steady_clock::now();

    // Watchdog: run on a separate thread, kill the test if it's still
    // grinding after 60 s.
    std::vector<Mas> results;
    auto fut = std::async(std::launch::async, [&]() {
        results = coilAdviser.get_advised_coil(mas, 1);
    });

    auto status = fut.wait_for(std::chrono::seconds(60));
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start).count();

    std::cout << "[HANG_REPRO] elapsed = " << elapsed << " ms" << std::endl;

    auto log = read_log();
    std::cout << "[HANG_REPRO] MKF log size = " << log.size() << std::endl;
    if (!log.empty()) {
        std::cout << "[HANG_REPRO] MKF log:\n" << log << std::endl;
    }

    if (status != std::future_status::ready) {
        std::cout << "[HANG_REPRO] WATCHDOG TIMEOUT — adviser did not return in 60 s"
                  << std::endl;
        // Detach the future; the std::async destructor would otherwise block
        // on the still-running task. We're about to abort the process via
        // REQUIRE anyway, so leaking the thread for the abort path is fine.
        std::quick_exit(124);  // 124 == GNU timeout sentinel
    }

    fut.get();
    std::cout << "[HANG_REPRO] results.size = " << results.size() << std::endl;
    REQUIRE(results.size() > 0);

    settings.reset();
}
