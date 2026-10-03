// Settings::allowMaterialDataExtrapolation and the Logger's structured collector.
//
// Off (the default) every use of a core material outside its data throws, as it always has. On, an
// explicit opt-in for debugging and prospecting, MKF evaluates the same model there and logs a WARNING
// naming the material, the quantity, the value and the range; the collector hands such records back to
// a caller (the web app shows them in its log panel).
#include "physical_models/CoreLosses.h"
#include "advisers/MagneticFilter.h"
#include "advisers/CoreAdviser.h"
#include "advisers/CoilAdviser.h"
#include "advisers/WireAdviser.h"
#include "advisers/MagneticAdviser.h"
#include "advisers/CoreCrossReferencer.h"
#include "advisers/CoreMaterialCrossReferencer.h"
#include "processors/Sweeper.h"
#include "support/Exceptions.h"
#include "support/Logger.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <numbers>
#include <thread>

using namespace OpenMagnetics;

namespace {

// Enables the process-wide collector for one test and always disables it again.
struct CollectorGuard {
    explicit CollectorGuard(LogLevel level = LogLevel::WARNING) {
        Logger::getInstance().disableCollector();
        Logger::getInstance().enableCollector(level);
    }
    ~CollectorGuard() { Logger::getInstance().disableCollector(); }
};

std::vector<LogRecord> extrapolation_records(const std::vector<LogRecord>& records) {
    std::vector<LogRecord> found;
    for (const auto& record : records) {
        if (record.moduleOfOrigin == kMaterialDataExtrapolationModule) {
            found.push_back(record);
        }
    }
    return found;
}

// One period of a sampled sinusoid (no time axis: the models scale it to the excitation's frequency),
// as TestCoreLosses' build_sinusoidal_flux_excitation does; the waveform-integrating models (iGSE, Roshen)
// need the samples, a processed-only sinusoid is not enough for them.
OperatingPointExcitation sinusoidal_flux_excitation(double frequency, double peak) {
    const size_t numberPoints = 10000;
    std::vector<double> data(numberPoints);
    for (size_t i = 0; i < numberPoints; ++i) {
        data[i] = peak * sin(2 * std::numbers::pi * static_cast<double>(i) / (numberPoints - 1));
    }
    json excitationJson;
    excitationJson["frequency"] = frequency;
    excitationJson["magneticFluxDensity"]["waveform"]["data"] = data;
    excitationJson["magneticFluxDensity"]["processed"]["label"] = WaveformLabel::SINUSOIDAL;
    excitationJson["magneticFluxDensity"]["processed"]["offset"] = 0;
    excitationJson["magneticFluxDensity"]["processed"]["peak"] = peak;
    excitationJson["magneticFluxDensity"]["processed"]["peakToPeak"] = 2 * peak;
    excitationJson["magneticFluxDensity"]["processed"]["dutyCycle"] = 0.5;
    return OperatingPointExcitation(excitationJson);
}

}  // namespace

TEST_CASE("Material data extrapolation is off by default and reset() turns it off", "[material-extrapolation][settings]") {
    settings.reset();
    CHECK_FALSE(settings.get_allow_material_data_extrapolation());
    settings.set_allow_material_data_extrapolation(true);
    CHECK(settings.get_allow_material_data_extrapolation());
    settings.reset();
    CHECK_FALSE(settings.get_allow_material_data_extrapolation());
}

TEST_CASE("Steinmetz outside the fitted span: off throws and logs nothing, on extrapolates and warns", "[material-extrapolation][core-losses]") {
    settings.reset();
    clear_databases();
    CollectorGuard collector;

    // TP5H is fitted over 1-5 MHz only; 100 kHz is below it.
    auto [tp5hMinimum, tp5hMaximum] = CoreLossesModel::get_steinmetz_fitted_span("TP5H");
    REQUIRE(tp5hMinimum > 100000);
    auto coreLossesModel = CoreLossesModel::factory(CoreLossesModels::STEINMETZ);
    auto excitation = sinusoidal_flux_excitation(100000, 0.05);

    SECTION("off") {
        CHECK_THROWS_AS(CoreLossesModel::get_steinmetz_coefficients("TP5H", 100000), MaterialFrequencyOutOfSpanException);
        CHECK_THROWS_AS(coreLossesModel->get_core_volumetric_losses(Core::resolve_material("TP5H"), excitation, 25),
                        MaterialFrequencyOutOfSpanException);
        CHECK(extrapolation_records(Logger::getInstance().drainCollected()).empty());
    }
    SECTION("on") {
        settings.set_allow_material_data_extrapolation(true);
        auto below = CoreLossesModel::get_steinmetz_coefficients("TP5H", 100000);
        // The range carried past its edge is the one nearest the frequency: the lowest.
        CHECK(below.get_minimum_frequency().value() == tp5hMinimum);
        double volumetricLosses = coreLossesModel->get_core_volumetric_losses(Core::resolve_material("TP5H"), excitation, 25);
        CHECK(std::isfinite(volumetricLosses));
        CHECK(volumetricLosses > 0);

        auto records = extrapolation_records(Logger::getInstance().drainCollected());
        REQUIRE_FALSE(records.empty());
        CHECK(records[0].level == LogLevel::WARNING);
        CHECK_THAT(records[0].message, Catch::Matchers::ContainsSubstring("TP5H") &&
                                       Catch::Matchers::ContainsSubstring("Steinmetz") &&
                                       Catch::Matchers::ContainsSubstring(std::to_string(100000.0)) &&
                                       Catch::Matchers::ContainsSubstring(std::to_string(tp5hMinimum)) &&
                                       Catch::Matchers::ContainsSubstring(std::to_string(tp5hMaximum)) &&
                                       Catch::Matchers::ContainsSubstring("below") &&
                                       Catch::Matchers::ContainsSubstring("allowMaterialDataExtrapolation"));

        // Above the highest range: 3C95's highest range is the one carried up.
        auto [c95Minimum, c95Maximum] = CoreLossesModel::get_steinmetz_fitted_span("3C95");
        auto above = CoreLossesModel::get_steinmetz_coefficients("3C95", 2 * c95Maximum);
        CHECK(above.get_maximum_frequency().value() == c95Maximum);
        auto aboveRecords = extrapolation_records(Logger::getInstance().drainCollected());
        REQUIRE(aboveRecords.size() == 1);
        CHECK_THAT(aboveRecords[0].message, Catch::Matchers::ContainsSubstring("3C95") && Catch::Matchers::ContainsSubstring("above"));

        // Inside the span nothing changes and nothing is logged.
        CHECK_NOTHROW(CoreLossesModel::get_steinmetz_coefficients("3C95", 100000));
        CHECK(extrapolation_records(Logger::getInstance().drainCollected()).empty());
    }
    settings.reset();
}

TEST_CASE("Steinmetz in a gap between ranges: off throws, on uses the nearest range and warns", "[material-extrapolation][core-losses]") {
    settings.reset();
    clear_databases();
    CollectorGuard collector;

    // 3C95 with its 150 kHz-1 MHz range moved to start at 200 kHz: 150-200 kHz is a gap.
    json materialJson;
    to_json(materialJson, Core::resolve_material("3C95"));
    bool gapMade = false;
    for (auto& method : materialJson["volumetricLosses"]["default"]) {
        if (method.is_object() && method.contains("method") && method["method"] == "steinmetz") {
            for (auto& range : method["ranges"]) {
                if (range["minimumFrequency"] == 150000.0) {
                    range["minimumFrequency"] = 200000.0;
                    gapMade = true;
                }
            }
        }
    }
    REQUIRE(gapMade);
    materialJson["name"] = "3C95 with a gap";
    CoreMaterial gapped;
    from_json(materialJson, gapped);

    SECTION("off") {
        CHECK_THROWS_WITH(CoreLossesModel::get_steinmetz_coefficients(gapped, 160000),
                          Catch::Matchers::ContainsSubstring("gap between its Steinmetz ranges"));
        CHECK(extrapolation_records(Logger::getInstance().drainCollected()).empty());
    }
    SECTION("on") {
        settings.set_allow_material_data_extrapolation(true);
        auto datum = CoreLossesModel::get_steinmetz_coefficients(gapped, 160000);
        CHECK(datum.get_maximum_frequency().value() == 150000.0);  // 10 kHz away, against 40 kHz to the next
        auto records = extrapolation_records(Logger::getInstance().drainCollected());
        REQUIRE(records.size() == 1);
        CHECK_THAT(records[0].message, Catch::Matchers::ContainsSubstring("3C95 with a gap") &&
                                       Catch::Matchers::ContainsSubstring("gap"));
    }
    settings.reset();
}

TEST_CASE("Core losses above the Curie temperature: off throws, on evaluates and warns", "[material-extrapolation][core-losses]") {
    settings.reset();
    clear_databases();
    CollectorGuard collector;
    Core core = OpenMagneticsTesting::get_quick_core("T 20/10/7", json::array(), 1, "PC95");
    double curieTemperature = core.resolve_material().get_curie_temperature().value();
    auto excitation = sinusoidal_flux_excitation(200000, 0.05);
    CoreLosses coreLosses;

    SECTION("off") {
        CHECK_THROWS_AS(coreLosses.get_core_volumetric_losses(core.resolve_material(), excitation, curieTemperature + 10),
                        MaterialAboveCurieTemperatureException);
        CHECK(extrapolation_records(Logger::getInstance().drainCollected()).empty());
    }
    SECTION("on") {
        settings.set_allow_material_data_extrapolation(true);
        double volumetricLosses = coreLosses.get_core_volumetric_losses(core.resolve_material(), excitation, curieTemperature + 10);
        CHECK(std::isfinite(volumetricLosses));
        CHECK(volumetricLosses > 0);
        auto records = extrapolation_records(Logger::getInstance().drainCollected());
        REQUIRE_FALSE(records.empty());
        CHECK(records[0].level == LogLevel::WARNING);
        CHECK_THAT(records[0].message, Catch::Matchers::ContainsSubstring("PC95") &&
                                       Catch::Matchers::ContainsSubstring("Curie") &&
                                       Catch::Matchers::ContainsSubstring(std::to_string(curieTemperature)) &&
                                       Catch::Matchers::ContainsSubstring(std::to_string(curieTemperature + 10)));
    }
    settings.reset();
}

TEST_CASE("Frequency from core losses sweeps outside the span only with the flag on", "[material-extrapolation][core-losses]") {
    settings.reset();
    clear_databases();
    CollectorGuard collector;
    // TP5H is fitted over 1-5 MHz; the loss to invert is the (extrapolated) iGSE loss at 500 kHz.
    Core core = OpenMagneticsTesting::get_quick_core("T 20/10/7", json::array(), 1, "TP5H");
    auto [spanMinimum, spanMaximum] = CoreLossesModel::get_steinmetz_fitted_span("TP5H");
    REQUIRE(spanMinimum > 500000);
    auto excitation = sinusoidal_flux_excitation(500000, 0.02);
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    auto coreLossesModel = CoreLossesModel::factory(CoreLossesModels::IGSE);

    settings.set_allow_material_data_extrapolation(true);
    double lossesAt500kHz = coreLossesModel->get_core_losses(core, excitation, 25).get_core_losses();
    REQUIRE(lossesAt500kHz > 0);
    Logger::getInstance().drainCollected();

    SECTION("on: the whole 10 kHz-2 MHz sweep runs and finds 500 kHz, with warnings") {
        double frequency = coreLossesModel->get_frequency_from_core_losses(core, magneticFluxDensity, 25, lossesAt500kHz);
        CHECK_THAT(frequency, Catch::Matchers::WithinAbs(500000, 5000));
        CHECK_FALSE(extrapolation_records(Logger::getInstance().drainCollected()).empty());
    }
    SECTION("off: the sweep is clipped to the span and logs nothing") {
        settings.set_allow_material_data_extrapolation(false);
        double frequency = coreLossesModel->get_frequency_from_core_losses(core, magneticFluxDensity, 25, lossesAt500kHz);
        CHECK(frequency >= spanMinimum);
        CHECK(extrapolation_records(Logger::getInstance().drainCollected()).empty());
    }
    settings.reset();
}

TEST_CASE("The loss-model span gate lets an out-of-span material through only with the flag on", "[material-extrapolation][adviser]") {
    settings.reset();
    auto inputs = OpenMagneticsTesting::create_quick_test_inputs();  // 100 kHz
    CHECK_FALSE(MagneticFilterLossModelFrequencySpan::is_material_evaluable(Core::resolve_material("TP5H"), inputs, CoreLossesModels::IGSE));
    settings.set_allow_material_data_extrapolation(true);
    CHECK(MagneticFilterLossModelFrequencySpan::is_material_evaluable(Core::resolve_material("TP5H"), inputs, CoreLossesModels::IGSE));
    // A material no model can evaluate stays not evaluable: the flag extrapolates data, it does not invent it.
    CHECK_FALSE(MagneticFilterLossModelFrequencySpan::is_material_evaluable(Core::resolve_material("NP7"), inputs));
    settings.reset();
}

// ABT #1652 (Alf): advisers NEVER extrapolate, whatever the flag says. Each adviser class holds a
// MaterialDataExtrapolationBarrier, so while one exists -- on this thread or a worker running on a Settings
// snapshot -- the flag reads false, the span gate bars the material and the coefficients throw. When the
// adviser is gone the caller's choice reads back unchanged.
template <typename Adviser>
static void check_adviser_bars_extrapolation(const std::string& adviserName, const OpenMagnetics::Inputs& inputs) {
    INFO(adviserName);
    REQUIRE(settings.get_allow_material_data_extrapolation());
    {
        Adviser adviser;
        CHECK_FALSE(settings.get_allow_material_data_extrapolation());
        CHECK_THROWS_AS(CoreLossesModel::get_steinmetz_coefficients("TP5H", 100000), MaterialFrequencyOutOfSpanException);
        CHECK_FALSE(MagneticFilterLossModelFrequencySpan::is_material_evaluable(Core::resolve_material("TP5H"), inputs, CoreLossesModels::IGSE));

        // A worker thread an adviser fans out to runs on a snapshot of the parent's Settings, which
        // carries the caller's flag: it must be barred as well.
        const Settings parentSnapshot = Settings::GetInstance();
        bool workerAllowed = true;
        bool workerThrew = false;
        std::thread worker([&]() {
            Settings::GetInstance() = parentSnapshot;
            workerAllowed = Settings::GetInstance().get_allow_material_data_extrapolation();
            try {
                CoreLossesModel::get_steinmetz_coefficients("TP5H", 100000);
            }
            catch (const MaterialFrequencyOutOfSpanException&) {
                workerThrew = true;
            }
        });
        worker.join();
        CHECK_FALSE(workerAllowed);
        CHECK(workerThrew);

        // A copy is an adviser too; destroying it does not lift the original's barrier.
        { Adviser copy(adviser); }
        CHECK_FALSE(settings.get_allow_material_data_extrapolation());
    }
    CHECK(settings.get_allow_material_data_extrapolation());
}

TEST_CASE("Advisers never extrapolate material data, whatever the flag says (ABT #1652)", "[material-extrapolation][adviser-barrier]") {
    settings.reset();
    clear_databases();
    CollectorGuard collector;
    auto inputs = OpenMagneticsTesting::create_quick_test_inputs();  // 100 kHz, below TP5H's 1-5 MHz fit
    REQUIRE(MaterialDataExtrapolationBarrier::active_count() == 0);
    settings.set_allow_material_data_extrapolation(true);
    REQUIRE(MagneticFilterLossModelFrequencySpan::is_material_evaluable(Core::resolve_material("TP5H"), inputs, CoreLossesModels::IGSE));
    Logger::getInstance().drainCollected();

    check_adviser_bars_extrapolation<CoreAdviser>("CoreAdviser", inputs);
    check_adviser_bars_extrapolation<WireAdviser>("WireAdviser", inputs);
    check_adviser_bars_extrapolation<CoilAdviser>("CoilAdviser", inputs);
    check_adviser_bars_extrapolation<MagneticAdviser>("MagneticAdviser", inputs);
    check_adviser_bars_extrapolation<CoreCrossReferencer>("CoreCrossReferencer", inputs);
    check_adviser_bars_extrapolation<CoreMaterialCrossReferencer>("CoreMaterialCrossReferencer", inputs);

    // Nested advisers (MagneticAdviser runs CoreAdviser runs CoilAdviser): the inner one ending does not
    // lift the outer one's barrier.
    {
        MagneticAdviser outer;
        { CoreAdviser inner; }
        CHECK_FALSE(settings.get_allow_material_data_extrapolation());
    }
    CHECK(MaterialDataExtrapolationBarrier::active_count() == 0);
    CHECK(settings.get_allow_material_data_extrapolation());
    // Nothing was extrapolated under any adviser.
    CHECK(extrapolation_records(Logger::getInstance().drainCollected()).empty());
    settings.reset();
}

TEST_CASE("The core-resistance sweep keeps the requested window only with the flag on", "[material-extrapolation][sweeper]") {
    settings.reset();
    clear_databases();
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic("E 42/21/15", json::array(), {10}, 1, "3C95");
    auto [c95Minimum, c95Maximum] = CoreLossesModel::get_steinmetz_fitted_span("3C95");
    auto clipped = Sweeper::core_resistance_frequency_window(magnetic, c95Minimum / 10, c95Maximum * 2);
    CHECK(clipped.first == c95Minimum);
    CHECK(clipped.second == c95Maximum);
    settings.set_allow_material_data_extrapolation(true);
    auto requested = Sweeper::core_resistance_frequency_window(magnetic, c95Minimum / 10, c95Maximum * 2);
    CHECK(requested.first == c95Minimum / 10);
    CHECK(requested.second == c95Maximum * 2);
    settings.reset();
}

TEST_CASE("Logger collector: own level, merged duplicates, drain empties, console level untouched", "[material-extrapolation][logger]") {
    auto& logger = Logger::getInstance();
    const LogLevel consoleLevel = logger.getLevel();
    logger.disableCollector();

    SECTION("disabled collects nothing") {
        logger.warning("not collected", "TestModule");
        logger.enableCollector(LogLevel::WARNING);
        CHECK(logger.drainCollected().empty());
        logger.disableCollector();
    }
    SECTION("level, merging and order") {
        CollectorGuard collector(LogLevel::INFO);
        CHECK(logger.isCollectorEnabled());
        CHECK(logger.getCollectorLevel() == LogLevel::INFO);
        logger.debug("below the collector level", "TestModule");
        logger.info("first", "TestModule");
        logger.warning("second", "OtherModule");
        logger.info("first", "TestModule");
        logger.error("third", "TestModule");
        auto records = logger.drainCollected();
        REQUIRE(records.size() == 3);
        CHECK(records[0].level == LogLevel::INFO);
        CHECK(records[0].moduleOfOrigin == "TestModule");
        CHECK(records[0].message == "first");
        CHECK(records[0].count == 2);
        CHECK(records[1].level == LogLevel::WARNING);
        CHECK(records[1].moduleOfOrigin == "OtherModule");
        CHECK(records[2].level == LogLevel::ERROR);
        CHECK(logger.drainCollected().empty());
        // The collector does not move the console sink's level.
        CHECK(logger.getLevel() == consoleLevel);
    }
    SECTION("cap") {
        CollectorGuard collector(LogLevel::WARNING);
        for (size_t index = 0; index < Logger::kMaximumCollectedRecords + 5; ++index) {
            logger.warning("record " + std::to_string(index), "TestModule");
        }
        auto records = logger.drainCollected();
        REQUIRE(records.size() == Logger::kMaximumCollectedRecords + 1);
        CHECK(records.back().moduleOfOrigin == "Logger");
        CHECK_THAT(records.back().message, Catch::Matchers::StartsWith("5 further log records were dropped"));
    }
    SECTION("disable discards") {
        logger.enableCollector(LogLevel::WARNING);
        logger.warning("discarded", "TestModule");
        logger.disableCollector();
        logger.enableCollector(LogLevel::WARNING);
        CHECK(logger.drainCollected().empty());
        logger.disableCollector();
    }
    SECTION("level names") {
        CHECK(log_level_from_string("WARNING") == LogLevel::WARNING);
        CHECK(log_level_from_string("OFF") == LogLevel::OFF);
        CHECK_THROWS_AS(log_level_from_string("warn"), std::invalid_argument);
    }
    CHECK(logger.getLevel() == consoleLevel);
}
