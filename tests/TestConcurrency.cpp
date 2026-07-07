// ABT #113: global-state reentrancy harness.
//
// MKF has historically been driven strictly single-threaded (see the
// THREAD-SAFETY note in support/Utils.h). This file is the validation
// instrument for making the compute paths reentrant: it fans the
// CoreAdviser (AVAILABLE_CORES, small fixed core shortlist) across several
// std::thread workers, each with its OWN Inputs, and asserts:
//   1. no crash / no exception escapes a worker,
//   2. every worker's result set is non-empty,
//   3. DETERMINISM: each worker's result is identical to the same query run
//      single-threaded (top-result core reference + score to 1e-9 rel).
//
// Tagged [concurrency][heavy] — NOT part of [smoke-test]. Catch2 assertion
// macros are not thread-safe, so workers record into pre-sized slots and the
// main thread does all REQUIREs after join.
#include "support/Settings.h"
#include "support/Utils.h"
#include "advisers/CoreAdviser.h"
#include "physical_models/WindingLosses.h"
#include "processors/Inputs.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>

#include <atomic>
#include <barrier>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

namespace {

constexpr size_t NUMBER_THREADS = 6;
constexpr size_t ITERATIONS_PER_THREAD = 3;
constexpr size_t CORE_SHORTLIST_SIZE = 60;
constexpr double SCORE_RELATIVE_TOLERANCE = 1e-9;

struct AdviserQuery {
    double desiredInductance;
    double frequency;
    double voltagePeakToPeak;
};

// Distinct queries so concurrent workers exercise different code paths /
// cache keys. Thread i runs query i % queries.size().
const std::vector<AdviserQuery> QUERIES = {
    {10e-5, 100000, 600},
    {5e-5, 200000, 400},
    {2e-5, 400000, 200},
    {20e-5, 50000, 800},
    {8e-5, 150000, 500},
    {15e-5, 75000, 700},
};

OpenMagnetics::Inputs make_inputs(const AdviserQuery& query) {
    return OpenMagnetics::Inputs::create_quick_operating_point(
        query.frequency, query.desiredInductance, 25, WaveformLabel::SINUSOIDAL,
        query.voltagePeakToPeak, 0.5, 0, {});
}

// Small fixed core shortlist for speed: first N parseable cores of the
// test_cores.ndjson fixture. Loaded once on the main thread; each worker
// gets its own COPY (get_advised_core moves elements out of the vector it
// is given).
std::vector<Core> load_core_shortlist() {
    auto inventoryPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_cores.ndjson");
    std::ifstream ndjsonFile(inventoryPath);
    std::string jsonLine;
    std::vector<Core> cores;
    while (cores.size() < CORE_SHORTLIST_SIZE && std::getline(ndjsonFile, jsonLine)) {
        json jf = json::parse(jsonLine);
        try {
            cores.emplace_back(jf, false, true, false);
        }
        catch (const CoreShapeNotFoundException&) {
            continue;
        }
    }
    return cores;
}

struct AdviserOutcome {
    bool completed = false;
    std::string error;
    size_t numberResults = 0;
    std::string topReference;
    double topScore = std::numeric_limits<double>::quiet_NaN();
};

AdviserOutcome run_core_adviser_query(const AdviserQuery& query, const std::vector<Core>& coreShortlist) {
    AdviserOutcome outcome;
    try {
        auto inputs = make_inputs(query);
        std::map<CoreAdviser::CoreAdviserFilters, double> weights;
        weights[CoreAdviser::CoreAdviserFilters::COST] = 1;
        weights[CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
        weights[CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        CoreAdviser coreAdviser;
        coreAdviser.set_mode(CoreAdviser::CoreAdviserModes::AVAILABLE_CORES);
        auto cores = coreShortlist;  // own mutable copy per call
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, &cores, 5);

        outcome.numberResults = masMagnetics.size();
        if (!masMagnetics.empty()) {
            outcome.topReference = masMagnetics[0].first.get_magnetic().get_core().get_name().value_or("<unnamed>");
            outcome.topScore = masMagnetics[0].second;
        }
        outcome.completed = true;
    }
    catch (const std::exception& e) {
        outcome.error = e.what();
    }
    catch (...) {
        outcome.error = "non-std exception";
    }
    return outcome;
}

bool scores_match(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) {
        return false;
    }
    double scale = std::max(std::abs(a), std::abs(b));
    return std::abs(a - b) <= SCORE_RELATIVE_TOLERANCE * std::max(scale, 1e-300);
}

} // namespace

TEST_CASE("Test_Concurrency_CoreAdviser_Parallel_Deterministic", "[concurrency][heavy]") {
    settings.reset();
    clear_databases();

    auto coreShortlist = load_core_shortlist();
    REQUIRE(coreShortlist.size() >= CORE_SHORTLIST_SIZE / 2);

    // The parallel region runs FIRST, on cold memo caches: concurrent workers
    // must build their interpolator/inductance memos (and, pre-fix, would race
    // to lazy-load the shared catalogs) instead of merely reading state warmed
    // by a preceding single-threaded pass. The single-threaded references are
    // computed AFTER the join; determinism therefore also proves the memo
    // caches are semantically transparent (cold vs warm gives the same
    // physics).
    // Pre-fix version of this harness: no freeze/force-load API exists yet,
    // so the workers themselves trigger the lazy catalog loads concurrently —
    // exactly the disease this harness must detect.

    std::vector<std::vector<AdviserOutcome>> outcomes(NUMBER_THREADS, std::vector<AdviserOutcome>(ITERATIONS_PER_THREAD));
    std::barrier startBarrier(NUMBER_THREADS);
    std::vector<std::thread> workers;
    for (size_t threadIndex = 0; threadIndex < NUMBER_THREADS; ++threadIndex) {
        workers.emplace_back([threadIndex, &outcomes, &startBarrier, &coreShortlist] {
            const auto& query = QUERIES[threadIndex % QUERIES.size()];
            for (size_t iteration = 0; iteration < ITERATIONS_PER_THREAD; ++iteration) {
                // Line every worker up before each iteration to maximize overlap.
                startBarrier.arrive_and_wait();
                outcomes[threadIndex][iteration] = run_core_adviser_query(query, coreShortlist);
            }
        });
    }
    for (auto& worker : workers) {
        worker.join();
    }

    // Single-threaded references, computed after the parallel region.
    std::vector<AdviserOutcome> references;
    for (auto& query : QUERIES) {
        auto reference = run_core_adviser_query(query, coreShortlist);
        INFO("single-threaded reference failed: " << reference.error);
        REQUIRE(reference.completed);
        REQUIRE(reference.numberResults > 0);
        references.push_back(reference);
    }
    // Repeat once single-threaded: the reference must be self-consistent,
    // otherwise the determinism assertion below would be meaningless.
    for (size_t queryIndex = 0; queryIndex < QUERIES.size(); ++queryIndex) {
        auto again = run_core_adviser_query(QUERIES[queryIndex], coreShortlist);
        REQUIRE(again.completed);
        REQUIRE(again.topReference == references[queryIndex].topReference);
        REQUIRE(scores_match(again.topScore, references[queryIndex].topScore));
    }

    for (size_t threadIndex = 0; threadIndex < NUMBER_THREADS; ++threadIndex) {
        const auto& reference = references[threadIndex % QUERIES.size()];
        for (size_t iteration = 0; iteration < ITERATIONS_PER_THREAD; ++iteration) {
            const auto& outcome = outcomes[threadIndex][iteration];
            INFO("thread " << threadIndex << " iteration " << iteration << " error: " << outcome.error);
            CHECK(outcome.completed);
            if (!outcome.completed) {
                continue;
            }
            CHECK(outcome.numberResults > 0);
            INFO("thread " << threadIndex << " iteration " << iteration
                 << " top=" << outcome.topReference << " score=" << outcome.topScore
                 << " expected top=" << reference.topReference << " score=" << reference.topScore);
            CHECK(outcome.topReference == reference.topReference);
            CHECK(scores_match(outcome.topScore, reference.topScore));
        }
    }
    settings.reset();
}

TEST_CASE("Test_Concurrency_CoreAdviser_With_WindingLosses_Mixed", "[concurrency][heavy]") {
    settings.reset();
    clear_databases();

    auto coreShortlist = load_core_shortlist();
    REQUIRE(coreShortlist.size() >= CORE_SHORTLIST_SIZE / 2);

    // Fixed magnetic (wound coil) + operating point for the WindingLosses thread.
    std::vector<OpenMagnetics::Wire> wires;
    {
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_diameter(0.00071);
        wire.set_nominal_value_outer_diameter(0.000762);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::ROUND);
        wires.push_back(wire);
    }
    auto coil = OpenMagneticsTesting::get_quick_coil({20}, {1}, "ETD 29/16/10", 1,
        WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
        CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires, true, 1);
    auto core = OpenMagneticsTesting::get_quick_core("ETD 29/16/10", OpenMagneticsTesting::get_ground_gap(0.001), 1, "N87");
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    auto lossesInputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(
        100000, 1e-3, 25, WaveformLabel::SINUSOIDAL, 2 * 1.4142, 0.5, 0);
    auto lossesOperatingPoint = lossesInputs.get_operating_point(0);

    double referenceLosses = WindingLosses().calculate_losses(magnetic, lossesOperatingPoint, 25).get_winding_losses();
    REQUIRE(referenceLosses > 0);

    constexpr size_t NUMBER_ADVISER_THREADS = 4;
    constexpr size_t LOSSES_ITERATIONS = 6;
    std::vector<AdviserOutcome> adviserReferences;
    for (size_t threadIndex = 0; threadIndex < NUMBER_ADVISER_THREADS; ++threadIndex) {
        auto reference = run_core_adviser_query(QUERIES[threadIndex % QUERIES.size()], coreShortlist);
        INFO("single-threaded reference failed: " << reference.error);
        REQUIRE(reference.completed);
        REQUIRE(reference.numberResults > 0);
        adviserReferences.push_back(reference);
    }

    std::vector<AdviserOutcome> adviserOutcomes(NUMBER_ADVISER_THREADS);
    std::vector<double> lossesResults(LOSSES_ITERATIONS, std::numeric_limits<double>::quiet_NaN());
    std::string lossesError;
    std::barrier startBarrier(NUMBER_ADVISER_THREADS + 1);

    std::vector<std::thread> workers;
    for (size_t threadIndex = 0; threadIndex < NUMBER_ADVISER_THREADS; ++threadIndex) {
        workers.emplace_back([threadIndex, &adviserOutcomes, &startBarrier, &coreShortlist] {
            startBarrier.arrive_and_wait();
            adviserOutcomes[threadIndex] = run_core_adviser_query(QUERIES[threadIndex % QUERIES.size()], coreShortlist);
        });
    }
    workers.emplace_back([&, LOSSES_ITERATIONS] {
        startBarrier.arrive_and_wait();
        try {
            for (size_t iteration = 0; iteration < LOSSES_ITERATIONS; ++iteration) {
                lossesResults[iteration] = WindingLosses().calculate_losses(magnetic, lossesOperatingPoint, 25).get_winding_losses();
            }
        }
        catch (const std::exception& e) {
            lossesError = e.what();
        }
        catch (...) {
            lossesError = "non-std exception";
        }
    });
    for (auto& worker : workers) {
        worker.join();
    }

    for (size_t threadIndex = 0; threadIndex < NUMBER_ADVISER_THREADS; ++threadIndex) {
        const auto& outcome = adviserOutcomes[threadIndex];
        const auto& reference = adviserReferences[threadIndex];
        INFO("adviser thread " << threadIndex << " error: " << outcome.error);
        CHECK(outcome.completed);
        if (!outcome.completed) {
            continue;
        }
        CHECK(outcome.numberResults > 0);
        CHECK(outcome.topReference == reference.topReference);
        CHECK(scores_match(outcome.topScore, reference.topScore));
    }
    INFO("winding-losses thread error: " << lossesError);
    CHECK(lossesError.empty());
    for (size_t iteration = 0; iteration < LOSSES_ITERATIONS; ++iteration) {
        INFO("winding losses iteration " << iteration << " -> " << lossesResults[iteration]
             << " expected " << referenceLosses);
        CHECK(scores_match(lossesResults[iteration], referenceLosses));
    }
    settings.reset();
}
