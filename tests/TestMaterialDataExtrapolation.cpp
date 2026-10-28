// Settings::allowMaterialDataExtrapolation and the Logger's structured collector.
//
// On (the default since Alf, 2026-10-25: "when the data comes from the MAS/user, we just extrapolate and
// throw a warning, not an error"), MKF evaluates a core material outside its data with the same model and
// logs a WARNING naming the material, the quantity, the value and the range; the collector hands such
// records back to a caller (the web app shows them in its log panel). Off, or under any adviser (ABT #1652:
// advisers never extrapolate), every such use throws.
#include "physical_models/CoreLosses.h"
#include "advisers/MagneticFilter.h"
#include "advisers/CoreAdviser.h"
#include "advisers/CoilAdviser.h"
#include "advisers/WireAdviser.h"
#include "advisers/MagneticAdviser.h"
#include "advisers/CoreCrossReferencer.h"
#include "advisers/CoreMaterialCrossReferencer.h"
#include "processors/Sweeper.h"
#include "processors/MagneticSimulator.h"
#include "support/Exceptions.h"
#include "support/Logger.h"
#include "support/Utils.h"
#include "support/NumberFormat.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <filesystem>
#include <source_location>
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

// Alf, 2026-10-25: user and MAS data extrapolate with a warning by default; advisers never (barrier).
TEST_CASE("Material data extrapolation is on by default and reset() turns it back on", "[material-extrapolation][settings]") {
    settings.reset();
    CHECK(settings.get_allow_material_data_extrapolation());
    settings.set_allow_material_data_extrapolation(false);
    CHECK_FALSE(settings.get_allow_material_data_extrapolation());
    settings.reset();
    CHECK(settings.get_allow_material_data_extrapolation());
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
        settings.set_allow_material_data_extrapolation(false);  // the default is on since Alf, 2026-10-25
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
        // The message names the material and the span its data covers; the frequency asked is the record's detail.
        CHECK_THAT(records[0].message, Catch::Matchers::ContainsSubstring("TP5H") &&
                                       Catch::Matchers::ContainsSubstring("Steinmetz") &&
                                       Catch::Matchers::ContainsSubstring(format_quantity(tp5hMinimum, "Hz")) &&
                                       Catch::Matchers::ContainsSubstring(format_quantity(tp5hMaximum, "Hz")) &&
                                       !Catch::Matchers::ContainsSubstring("allowMaterialDataExtrapolation"));
        CHECK_THAT(records[0].detail, Catch::Matchers::ContainsSubstring(format_quantity(100000.0, "Hz")) &&
                                      Catch::Matchers::ContainsSubstring("below"));

        // Above the highest range: 3C95's highest range is the one carried up.
        auto [c95Minimum, c95Maximum] = CoreLossesModel::get_steinmetz_fitted_span("3C95");
        auto above = CoreLossesModel::get_steinmetz_coefficients("3C95", 2 * c95Maximum);
        CHECK(above.get_maximum_frequency().value() == c95Maximum);
        auto aboveRecords = extrapolation_records(Logger::getInstance().drainCollected());
        REQUIRE(aboveRecords.size() == 1);
        CHECK_THAT(aboveRecords[0].message, Catch::Matchers::ContainsSubstring("3C95"));
        CHECK_THAT(aboveRecords[0].detail, Catch::Matchers::ContainsSubstring("above"));

        // Inside the span nothing changes and nothing is logged.
        CHECK_NOTHROW(CoreLossesModel::get_steinmetz_coefficients("3C95", 100000));
        CHECK(extrapolation_records(Logger::getInstance().drainCollected()).empty());
    }
    settings.reset();
}

// Repeats of one material's extrapolation at different values are one record: the message is the same for
// every value asked (material and the span its data covers), the value asked is the detail, so the web's
// Engine log shows one row with a count instead of one row per frequency.
TEST_CASE("Extrapolation warnings for one material merge into one record whatever value is asked", "[material-extrapolation][core-losses][logger]") {
    settings.reset();
    clear_databases();
    CollectorGuard collector;
    settings.set_allow_material_data_extrapolation(true);
    auto [tp5hMinimum, tp5hMaximum] = CoreLossesModel::get_steinmetz_fitted_span("TP5H");
    REQUIRE(tp5hMinimum > 200000);
    CoreLossesModel::get_steinmetz_coefficients("TP5H", 100000);
    CoreLossesModel::get_steinmetz_coefficients("TP5H", 200000);
    auto records = extrapolation_records(Logger::getInstance().drainCollected());
    REQUIRE(records.size() == 1);
    CHECK(records[0].count == 2);
    CHECK_THAT(records[0].detail, Catch::Matchers::ContainsSubstring(format_quantity(200000.0, "Hz")));
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
        settings.set_allow_material_data_extrapolation(false);  // the default is on since Alf, 2026-10-25
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
        settings.set_allow_material_data_extrapolation(false);  // the default is on since Alf, 2026-10-25
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
                                       Catch::Matchers::ContainsSubstring(format_significant(curieTemperature) + " C"));
        CHECK_THAT(records[0].detail, Catch::Matchers::ContainsSubstring(format_significant(curieTemperature + 10) + " C"));
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
    settings.set_allow_material_data_extrapolation(false);  // the default is on since Alf, 2026-10-25
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

// Settings::simulatorLossTemperatureLoop is for simulate and the builder: every adviser class, its copies and the
// worker threads it fans out to (on a snapshot of the caller's Settings) read it ON while the adviser exists, through
// the same barrier. The adviser gone, the caller's choice reads back.
template <typename Adviser>
static void check_adviser_keeps_loss_temperature_loop(const std::string& adviserName) {
    INFO(adviserName);
    REQUIRE_FALSE(settings.get_simulator_loss_temperature_loop());
    {
        Adviser adviser;
        CHECK(settings.get_simulator_loss_temperature_loop());
        const Settings parentSnapshot = Settings::GetInstance();
        bool workerLoop = false;
        std::thread worker([&]() {
            Settings::GetInstance() = parentSnapshot;
            workerLoop = Settings::GetInstance().get_simulator_loss_temperature_loop();
        });
        worker.join();
        CHECK(workerLoop);
        { Adviser copy(adviser); }
        CHECK(settings.get_simulator_loss_temperature_loop());
    }
    CHECK_FALSE(settings.get_simulator_loss_temperature_loop());
}

TEST_CASE("Advisers keep the loss-temperature loop, whatever simulatorLossTemperatureLoop says", "[adviser-barrier][abt-2003]") {
    settings.reset();
    REQUIRE(MaterialDataExtrapolationBarrier::active_count() == 0);
    struct Restore {
        ~Restore() { Settings::GetInstance().reset(); }
    } restore;
    settings.set_simulator_loss_temperature_loop(false);
    check_adviser_keeps_loss_temperature_loop<CoreAdviser>("CoreAdviser");
    check_adviser_keeps_loss_temperature_loop<WireAdviser>("WireAdviser");
    check_adviser_keeps_loss_temperature_loop<CoilAdviser>("CoilAdviser");
    check_adviser_keeps_loss_temperature_loop<MagneticAdviser>("MagneticAdviser");
    check_adviser_keeps_loss_temperature_loop<CoreCrossReferencer>("CoreCrossReferencer");
    check_adviser_keeps_loss_temperature_loop<CoreMaterialCrossReferencer>("CoreMaterialCrossReferencer");
}

TEST_CASE("The core-resistance sweep keeps the requested window only with the flag on", "[material-extrapolation][sweeper]") {
    settings.reset();
    clear_databases();
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic("E 42/21/15", json::array(), {10}, 1, "3C95");
    auto [c95Minimum, c95Maximum] = CoreLossesModel::get_steinmetz_fitted_span("3C95");
    settings.set_allow_material_data_extrapolation(false);  // the default is on since Alf, 2026-10-25
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
        logger.warning("not collected", "Coil");
        logger.enableCollector(LogLevel::WARNING);
        CHECK(logger.drainCollected().empty());
        logger.disableCollector();
    }
    SECTION("level, merging and order") {
        CollectorGuard collector(LogLevel::INFO);
        CHECK(logger.isCollectorEnabled());
        CHECK(logger.getCollectorLevel() == LogLevel::INFO);
        logger.debug("below the collector level", "Coil");
        logger.info("first", "Coil");
        logger.warning("second", "Inputs");
        logger.info("first", "Coil");
        logger.error("third", "Coil");
        auto records = logger.drainCollected();
        REQUIRE(records.size() == 3);
        CHECK(records[0].level == LogLevel::INFO);
        CHECK(records[0].moduleOfOrigin == "Winding");
        CHECK(records[0].message == "first");
        CHECK(records[0].count == 2);
        CHECK(records[1].level == LogLevel::WARNING);
        CHECK(records[1].moduleOfOrigin == "Operating points");
        CHECK(records[2].level == LogLevel::ERROR);
        CHECK(logger.drainCollected().empty());
        // The collector does not move the console sink's level.
        CHECK(logger.getLevel() == consoleLevel);
    }
    SECTION("cap") {
        CollectorGuard collector(LogLevel::WARNING);
        for (size_t index = 0; index < Logger::kMaximumCollectedRecords + 5; ++index) {
            logger.warning("record " + std::to_string(index), "Coil");
        }
        auto records = logger.drainCollected();
        REQUIRE(records.size() == Logger::kMaximumCollectedRecords + 1);
        CHECK(records.back().moduleOfOrigin == "Engine log");
        CHECK_THAT(records.back().message, Catch::Matchers::StartsWith("5 further log records were dropped"));
    }
    SECTION("disable discards") {
        logger.enableCollector(LogLevel::WARNING);
        logger.warning("discarded", "Coil");
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

// Alf, 2026-10-25: "01_1Hz: when the data comes from the MAS/user, we just extrapolate and throw a warning, not an
// error." With DEFAULT settings (no flag touched), a user's MAS whose operating point lies below its material's
// fitted Steinmetz span simulates to finite core losses and reports the extrapolation once, naming the material,
// the frequency and the span. Under an adviser (ABT #1652) the same MAS still throws.
TEST_CASE("A MAS below its material's fitted span extrapolates by default with one warning; advisers still throw", "[material-extrapolation][core-losses][default-extrapolation]") {
    settings.reset();
    clear_databases();
    REQUIRE(MaterialDataExtrapolationBarrier::active_count() == 0);
    auto examplePath = std::filesystem::path{std::source_location::current().file_name()}.parent_path()
                           .append("..").append("MAS").append("examples").append("01_simple_inductor_etd34_n87_1Hz.json");
    auto mas = OpenMagneticsTesting::mas_loader(examplePath.string());
    REQUIRE(mas.get_inputs().get_operating_points().size() == 1);
    const double frequency = mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_frequency();
    REQUIRE(frequency == 1.0);

    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    const std::string materialName = magnetic.get_core().get_material_name();
    REQUIRE(materialName == "N87");
    auto [spanMinimum, spanMaximum] = CoreLossesModel::get_steinmetz_fitted_span(magnetic.get_core().resolve_material());
    REQUIRE(frequency < spanMinimum);

    {
        CollectorGuard collector;
        OpenMagnetics::Mas simulated;
        REQUIRE_NOTHROW(simulated = MagneticSimulator().simulate(mas.get_inputs(), magnetic, true));
        REQUIRE(simulated.get_outputs().size() == 1);
        REQUIRE(simulated.get_outputs()[0].get_core_losses());
        double coreLosses = simulated.get_outputs()[0].get_core_losses()->get_core_losses();
        CHECK(std::isfinite(coreLosses));
        CHECK(coreLosses > 0);

        auto records = extrapolation_records(Logger::getInstance().drainCollected());
        REQUIRE(records.size() == 1);
        CHECK(records[0].level == LogLevel::WARNING);
        CHECK_THAT(records[0].message, Catch::Matchers::ContainsSubstring(materialName) &&
                                       Catch::Matchers::ContainsSubstring(format_quantity(spanMinimum, "Hz")) &&
                                       Catch::Matchers::ContainsSubstring(format_quantity(spanMaximum, "Hz")));
        CHECK_THAT(records[0].detail, Catch::Matchers::ContainsSubstring("asked at " + format_quantity(frequency, "Hz") + ", below"));
    }
    {
        CollectorGuard collector;
        CoreAdviser adviser;
        CHECK_FALSE(settings.get_allow_material_data_extrapolation());
        CHECK_THROWS_AS(MagneticSimulator().simulate(mas.get_inputs(), magnetic, true), MaterialFrequencyOutOfSpanException);
        CHECK(extrapolation_records(Logger::getInstance().drainCollected()).empty());
    }
    CHECK(settings.get_allow_material_data_extrapolation());
    settings.reset();
}
