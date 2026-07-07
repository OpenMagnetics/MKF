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
#include "support/LibraryContext.h"
#include "advisers/CoreAdviser.h"
#include "advisers/MagneticAdviser.h"
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

// RAII: guarantee the process is unfrozen on every exit path (a failing
// assertion between freeze and unfreeze must not poison later test cases,
// whose Settings listener calls reset() -> clear_loaded_cores()).
struct DatabasesFreezeGuard {
    DatabasesFreezeGuard() { set_databases_frozen(true); }
    ~DatabasesFreezeGuard() { set_databases_frozen(false); }
    DatabasesFreezeGuard(const DatabasesFreezeGuard&) = delete;
    DatabasesFreezeGuard& operator=(const DatabasesFreezeGuard&) = delete;
};

// Copy the spawning thread's Settings into this (worker) thread.
// Settings is a thread_local singleton: worker threads start with a
// default-constructed instance and must inherit the parent's configuration
// explicitly (see Settings::GetInstance docs).
void inherit_settings(const Settings& parentSnapshot) {
    Settings::GetInstance() = parentSnapshot;
}

// ABT #164: a single-family shape constraint (allowlist of exactly one family).
AdviserConstraints shape_family_constraint(const std::string& family) {
    AdviserConstraints constraints;
    constraints.shapeFamily.allowed = {family};
    return constraints;
}

struct MagneticAdviserOutcome {
    bool completed = false;
    std::string error;
    size_t numberResults = 0;
    bool allResultsRespectConstraint = false;  // every returned core is in the requested family
    std::string topReference;
    double topScore = std::numeric_limits<double>::quiet_NaN();
};

// Run one full MagneticAdviser design (STANDARD_CORES, POWER) under a
// per-call shape-family constraint. Each call constructs its OWN MagneticAdviser
// so the per-call `_constraints` is thread-local by construction; the fix under
// test threads that constraint to the nested CoreAdviser/CoilAdviser as LOCAL
// filtered views instead of swapping the process-shared catalogs.
MagneticAdviserOutcome run_magnetic_adviser_query(const AdviserQuery& query,
                                                  const AdviserConstraints& constraints) {
    MagneticAdviserOutcome outcome;
    try {
        auto inputs = make_inputs(query);
        MagneticAdviser adviser;
        adviser.set_application(MAS::MagneticApplication::POWER);
        adviser.set_core_mode(CoreAdviser::CoreAdviserModes::STANDARD_CORES);
        auto results = adviser.get_advised_magnetic(inputs, size_t(1), nullptr, constraints);

        outcome.numberResults = results.size();
        outcome.allResultsRespectConstraint = true;
        for (auto& [mas, score] : results) {
            auto family = mas.get_mutable_magnetic().get_core().get_shape_family();
            if (!constraints.shapeFamily.empty()
                && !acceptsCoreShapeFamily(constraints.shapeFamily, family)) {
                outcome.allResultsRespectConstraint = false;
            }
        }
        if (!results.empty()) {
            outcome.topReference = results[0].first.get_mutable_magnetic().get_core().get_name().value_or("<unnamed>");
            outcome.topScore = results[0].second;
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

} // namespace

TEST_CASE("Test_Concurrency_CoreAdviser_Parallel_Deterministic", "[concurrency][heavy]") {
    settings.reset();
    clear_databases();

    auto coreShortlist = load_core_shortlist();
    REQUIRE(coreShortlist.size() >= CORE_SHORTLIST_SIZE / 2);

    // The parallel region runs FIRST, on cold per-thread memo caches:
    // concurrent workers must build their interpolator/inductance memos
    // instead of merely reading state warmed by a preceding single-threaded
    // pass (the pre-fix version of this harness ran exactly like this and
    // SIGSEGVed on the then-shared caches/lazy catalog loads). The
    // single-threaded references are computed AFTER the join; determinism
    // therefore also proves the memo caches are semantically transparent
    // (cold vs warm gives the same physics).
    //
    // Read-only catalogs are force-loaded before fan-out and frozen for the
    // duration: any mutation attempt inside the parallel region throws.
    load_all_databases();
    const Settings parentSnapshot = Settings::GetInstance();
    DatabasesFreezeGuard freezeGuard;

    std::vector<std::vector<AdviserOutcome>> outcomes(NUMBER_THREADS, std::vector<AdviserOutcome>(ITERATIONS_PER_THREAD));
    std::barrier startBarrier(NUMBER_THREADS);
    std::vector<std::thread> workers;
    for (size_t threadIndex = 0; threadIndex < NUMBER_THREADS; ++threadIndex) {
        workers.emplace_back([threadIndex, &outcomes, &startBarrier, &coreShortlist, &parentSnapshot] {
            inherit_settings(parentSnapshot);
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
    set_databases_frozen(false);  // references below may lazy-load/clear freely

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

    load_all_databases();

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

    const Settings parentSnapshot = Settings::GetInstance();
    DatabasesFreezeGuard freezeGuard;

    std::vector<AdviserOutcome> adviserOutcomes(NUMBER_ADVISER_THREADS);
    std::vector<double> lossesResults(LOSSES_ITERATIONS, std::numeric_limits<double>::quiet_NaN());
    std::string lossesError;
    std::barrier startBarrier(NUMBER_ADVISER_THREADS + 1);

    std::vector<std::thread> workers;
    for (size_t threadIndex = 0; threadIndex < NUMBER_ADVISER_THREADS; ++threadIndex) {
        workers.emplace_back([threadIndex, &adviserOutcomes, &startBarrier, &coreShortlist, &parentSnapshot] {
            inherit_settings(parentSnapshot);
            startBarrier.arrive_and_wait();
            adviserOutcomes[threadIndex] = run_core_adviser_query(QUERIES[threadIndex % QUERIES.size()], coreShortlist);
        });
    }
    workers.emplace_back([&, LOSSES_ITERATIONS] {
        inherit_settings(parentSnapshot);
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

// ABT #164 (1): regression for the CoreAdviser batching path gutting the
// caller's vector. The batching overload used to std::move elements OUT of the
// vector it was handed; get_advised_core(...) feeds it a pointer to the
// process-SHARED coreDatabase, so a single batched run emptied the global core
// catalog for every subsequent caller (a latent single-threaded bug, and a
// silent corruption of the frozen fan-out catalog). Single-threaded on purpose;
// not [heavy].
TEST_CASE("Test_Concurrency_CoreAdviser_Batching_Preserves_Shared_Database", "[concurrency][core-adviser]") {
    settings.reset();
    clear_databases();

    // Populate the SHARED catalog with a small, fast shortlist so the batched
    // run is quick while still exercising the exact hazard: batching over
    // &coreDatabase.
    coreDatabase = load_core_shortlist();
    const size_t originalSize = coreDatabase.size();
    REQUIRE(originalSize >= CORE_SHORTLIST_SIZE / 2);

    auto inputs = make_inputs(QUERIES[0]);
    std::map<CoreAdviser::CoreAdviserFilters, double> weights;
    weights[CoreAdviser::CoreAdviserFilters::COST] = 1;
    weights[CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
    weights[CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

    CoreAdviser coreAdviser;
    coreAdviser.set_mode(CoreAdviser::CoreAdviserModes::AVAILABLE_CORES);

    // Small batch size => several move-into-batches iterations against the
    // shared vector (the pre-fix code emptied it here).
    const size_t maximumNumberCores = std::max<size_t>(1, originalSize / 4);
    auto first = coreAdviser.get_advised_core(inputs, weights, &coreDatabase, size_t(5), maximumNumberCores);

    // The shared catalog must be intact after a batched run.
    CHECK(coreDatabase.size() == originalSize);
    REQUIRE(!first.empty());
    const std::string firstTop = first[0].first.get_magnetic().get_core().get_name().value_or("<unnamed>");
    const double firstScore = first[0].second;

    // A second identical run must return identical top results — only possible
    // if the shared database survived the first run intact.
    auto second = coreAdviser.get_advised_core(inputs, weights, &coreDatabase, size_t(5), maximumNumberCores);
    CHECK(coreDatabase.size() == originalSize);
    REQUIRE(!second.empty());
    CHECK(second.size() == first.size());
    CHECK(second[0].first.get_magnetic().get_core().get_name().value_or("<unnamed>") == firstTop);
    CHECK(scores_match(second[0].second, firstScore));

    clear_databases();  // drop the injected shortlist so later tests reload the real catalog
    settings.reset();
}

// ABT #164 (2): MagneticAdviser fan-out. Mirrors the CoreAdviser parallel test
// but exercises the FULL adviser (CoreAdviser -> CoilAdviser -> analytical
// MagneticSimulator -> scoring) with a per-thread shape-family constraint.
// Before the fix, the ctx-aware overloads swapped the process-SHARED
// coreDatabase/wireDatabase for a constraint-filtered copy for the duration of
// the call, so concurrent workers saw each other's filtered catalogs. With the
// fix each worker threads its constraint as a LOCAL view; the determinism check
// (concurrent == isolated single-threaded reference) plus the per-worker
// family check prove the isolation.
TEST_CASE("Test_Concurrency_MagneticAdviser_Parallel_Constraints_Isolated", "[concurrency][heavy]") {
    settings.reset();
    clear_databases();

    // Distinct (shape-family, query) specs per worker. ETD and PQ select
    // disjoint cores, so a leaked/shared filtered catalog would surface as
    // either a wrong-family core or a mismatch against the isolated
    // single-threaded reference. Each family is paired with a query known to
    // yield a feasible design under that family so the fan-out never degenerates
    // to empty results (the test asserts >0 to prove the workers actually ran).
    struct ThreadSpec { std::string family; size_t queryIndex; };
    const std::vector<ThreadSpec> THREAD_SPECS = {
        {"ETD", 0},
        {"PQ",  1},
        {"ETD", 0},
        {"PQ",  1},
    };
    const size_t NUMBER_MAGNETIC_THREADS = THREAD_SPECS.size();
    constexpr size_t MAGNETIC_ITERATIONS_PER_THREAD = 2;

    load_all_databases();
    const Settings parentSnapshot = Settings::GetInstance();

    // Single-threaded references (isolated), computed BEFORE the parallel region
    // so the fan-out determinism check compares against a known-good baseline.
    std::vector<MagneticAdviserOutcome> references(NUMBER_MAGNETIC_THREADS);
    for (size_t threadIndex = 0; threadIndex < NUMBER_MAGNETIC_THREADS; ++threadIndex) {
        const auto& spec = THREAD_SPECS[threadIndex];
        references[threadIndex] = run_magnetic_adviser_query(QUERIES[spec.queryIndex], shape_family_constraint(spec.family));
        INFO("single-threaded MagneticAdviser reference " << threadIndex << " failed: " << references[threadIndex].error);
        REQUIRE(references[threadIndex].completed);
        REQUIRE(references[threadIndex].numberResults > 0);
        REQUIRE(references[threadIndex].allResultsRespectConstraint);
    }
    // Sanity: the ETD-constrained and PQ-constrained references must select
    // different cores — otherwise the constraint isn't biting and the isolation
    // assertions below would be vacuous.
    REQUIRE(references[0].topReference != references[1].topReference);

    // Freeze the read-only catalogs for the parallel region; ANY shared-catalog
    // mutation inside a worker (the pre-fix swap did exactly this) is then a
    // loud failure or a determinism mismatch.
    DatabasesFreezeGuard freezeGuard;

    std::vector<std::vector<MagneticAdviserOutcome>> outcomes(
        NUMBER_MAGNETIC_THREADS, std::vector<MagneticAdviserOutcome>(MAGNETIC_ITERATIONS_PER_THREAD));
    std::barrier startBarrier(NUMBER_MAGNETIC_THREADS);
    std::vector<std::thread> workers;
    for (size_t threadIndex = 0; threadIndex < NUMBER_MAGNETIC_THREADS; ++threadIndex) {
        workers.emplace_back([threadIndex, &outcomes, &startBarrier, &parentSnapshot, &THREAD_SPECS] {
            inherit_settings(parentSnapshot);
            const auto& spec = THREAD_SPECS[threadIndex];
            const auto& query = QUERIES[spec.queryIndex];
            const auto constraint = shape_family_constraint(spec.family);
            for (size_t iteration = 0; iteration < MAGNETIC_ITERATIONS_PER_THREAD; ++iteration) {
                startBarrier.arrive_and_wait();
                outcomes[threadIndex][iteration] = run_magnetic_adviser_query(query, constraint);
            }
        });
    }
    for (auto& worker : workers) {
        worker.join();
    }
    set_databases_frozen(false);

    for (size_t threadIndex = 0; threadIndex < NUMBER_MAGNETIC_THREADS; ++threadIndex) {
        const auto& reference = references[threadIndex];
        for (size_t iteration = 0; iteration < MAGNETIC_ITERATIONS_PER_THREAD; ++iteration) {
            const auto& outcome = outcomes[threadIndex][iteration];
            INFO("MagneticAdviser thread " << threadIndex << " iteration " << iteration
                 << " (family=" << THREAD_SPECS[threadIndex].family << ") error: " << outcome.error);
            CHECK(outcome.completed);
            if (!outcome.completed) {
                continue;
            }
            CHECK(outcome.numberResults > 0);
            CHECK(outcome.allResultsRespectConstraint);
            INFO("thread " << threadIndex << " iteration " << iteration
                 << " top=" << outcome.topReference << " score=" << outcome.topScore
                 << " expected top=" << reference.topReference << " score=" << reference.topScore);
            CHECK(outcome.topReference == reference.topReference);
            CHECK(scores_match(outcome.topScore, reference.topScore));
        }
    }
    settings.reset();
}
