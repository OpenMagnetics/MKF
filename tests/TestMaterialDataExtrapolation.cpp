// Settings::allowMaterialDataExtrapolation and the Logger's structured collector.
//
// On (the default since Alf, 2026-10-25: "when the data comes from the MAS/user, we just extrapolate and
// throw a warning, not an error"), MKF evaluates a core material outside its data with the same model and
// logs a WARNING naming the material, the quantity, the value and the range; the collector hands such
// records back to a caller (the web app shows them in its log panel). Off, or under any adviser (ABT #1652:
// advisers never extrapolate), every such use throws.
#include "physical_models/CoreLosses.h"
#include "physical_models/ComplexPermeability.h"
#include "physical_models/InitialPermeability.h"
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
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <map>
#include <set>
#include <filesystem>
#include <functional>
#include <optional>
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

// ABT #2235 (Alf, 2026-11-02: "make sure that all material admit extrapolation when the proper flag is set, not
// just throw"). Complex permeability of a material the estimator refuses (here an iron powder): outside its table
// there is no estimate, so off it throws; on, each part is held at its nearest tabulated value, with a WARNING.
TEST_CASE("Complex permeability the estimator refuses: off throws, on holds the nearest tabulated value and warns", "[material-extrapolation][complex-permeability]") {
    settings.reset();
    clear_databases();
    CollectorGuard collector;
    auto material = Core::resolve_material("Mix 26");
    REQUIRE(material.get_permeability().get_complex());
    ComplexPermeability complexPermeability;
    auto [minimumFrequency, maximumFrequency] = complexPermeability.get_frequency_range(material);
    REQUIRE(minimumFrequency > 100);
    // The premise: the estimator has no estimate for this material (a refusal, not an estimate).
    REQUIRE(complexPermeability.get_estimated_ranges(material).empty());

    SECTION("off") {
        settings.set_allow_material_data_extrapolation(false);
        CHECK_THROWS_AS(complexPermeability.get_complex_permeability(material, minimumFrequency / 10),
                        ComplexPermeabilityEstimateUnavailableException);
        CHECK_THROWS_AS(complexPermeability.get_complex_permeability(material, maximumFrequency * 10),
                        ComplexPermeabilityEstimateUnavailableException);
        CHECK(extrapolation_records(Logger::getInstance().drainCollected()).empty());
    }
    SECTION("on") {
        settings.set_allow_material_data_extrapolation(true);
        auto lowEdge = complexPermeability.get_tabulated_complex_permeability(material, minimumFrequency);
        auto below = complexPermeability.get_complex_permeability_with_source(material, minimumFrequency / 10);
        CHECK(below.real == lowEdge.first);
        CHECK(below.imaginary == lowEdge.second);
        CHECK(below.realSource == ComplexPermeabilitySource::HELD_FLAT);
        CHECK(below.imaginarySource == ComplexPermeabilitySource::HELD_FLAT);
        CHECK(below.is_estimated());
        auto records = extrapolation_records(Logger::getInstance().drainCollected());
        REQUIRE_FALSE(records.empty());
        CHECK(records[0].level == LogLevel::WARNING);
        CHECK_THAT(records[0].message, Catch::Matchers::ContainsSubstring("Mix 26") &&
                                       Catch::Matchers::ContainsSubstring("complex permeability") &&
                                       Catch::Matchers::ContainsSubstring(format_quantity(minimumFrequency, "Hz")));
        CHECK_THAT(records[0].detail, Catch::Matchers::ContainsSubstring("below"));
        // Above: each part is held at the end of its own table (one part may be tabulated further than the
        // range where both are), so the value no longer changes with frequency.
        auto above = complexPermeability.get_complex_permeability_with_source(material, maximumFrequency * 10);
        auto farAbove = complexPermeability.get_complex_permeability_with_source(material, maximumFrequency * 1000);
        CHECK(above.realSource == ComplexPermeabilitySource::HELD_FLAT);
        CHECK(above.imaginarySource == ComplexPermeabilitySource::HELD_FLAT);
        CHECK(farAbove.real == above.real);
        CHECK(farAbove.imaginary == above.imaginary);
        CHECK(std::isfinite(above.real));
        CHECK(std::isfinite(above.imaginary));
        auto aboveRecords = extrapolation_records(Logger::getInstance().drainCollected());
        REQUIRE_FALSE(aboveRecords.empty());
        CHECK_THAT(aboveRecords[0].detail, Catch::Matchers::ContainsSubstring("above"));
        // Inside the table nothing is logged.
        complexPermeability.get_complex_permeability(material, std::sqrt(minimumFrequency * maximumFrequency));
        CHECK(extrapolation_records(Logger::getInstance().drainCollected()).empty());
    }
    settings.reset();
}

// ABT #2235: remanence, coercive force and saturation field strength run to zero at the Curie temperature past
// their last data point. A material with no Curie temperature has no such end: off it throws, on the property is
// held at its last tabulated value, with a WARNING. A material without the property at all still throws.
TEST_CASE("Ferromagnetic properties above their data without a Curie temperature: off throws, on holds the last value and warns", "[material-extrapolation][saturation]") {
    settings.reset();
    clear_databases();
    CollectorGuard collector;
    json materialJson;
    to_json(materialJson, Core::resolve_material("3C95"));
    materialJson.erase("curieTemperature");
    materialJson["name"] = "3C95 without a Curie temperature";
    CoreMaterial material;
    from_json(materialJson, material);
    REQUIRE_FALSE(material.get_curie_temperature());
    REQUIRE(material.get_remanence());
    double lastRemanenceTemperature = std::numeric_limits<double>::lowest();
    double lastRemanence = 0;
    auto remanenceData = material.get_remanence().value();  // a local: range-for over .value() of a temporary dangles
    for (const auto& datum : remanenceData) {
        if (datum.get_temperature() > lastRemanenceTemperature) {
            lastRemanenceTemperature = datum.get_temperature();
            lastRemanence = datum.get_magnetic_flux_density();
        }
    }
    // Above every tabulated point of the three properties (3C95's data ends at 100 C).
    double temperature = 300;
    REQUIRE(temperature > lastRemanenceTemperature);

    SECTION("off") {
        settings.set_allow_material_data_extrapolation(false);
        CHECK_THROWS_AS(Core::get_remanence(material, temperature), MaterialException);
        CHECK_THROWS_AS(Core::get_coercive_force(material, temperature), MaterialException);
        CHECK_THROWS_AS(Core::get_magnetic_field_strength_saturation(material, temperature), MaterialException);
        CHECK(extrapolation_records(Logger::getInstance().drainCollected()).empty());
    }
    SECTION("on") {
        settings.set_allow_material_data_extrapolation(true);
        CHECK(Core::get_remanence(material, temperature) == lastRemanence);
        CHECK(std::isfinite(Core::get_coercive_force(material, temperature)));
        CHECK(std::isfinite(Core::get_magnetic_field_strength_saturation(material, temperature)));
        auto records = extrapolation_records(Logger::getInstance().drainCollected());
        REQUIRE(records.size() == 3);
        CHECK_THAT(records[0].message, Catch::Matchers::ContainsSubstring("3C95 without a Curie temperature") &&
                                       Catch::Matchers::ContainsSubstring("remanence") &&
                                       Catch::Matchers::ContainsSubstring("no Curie temperature"));
        CHECK_THAT(records[0].detail, Catch::Matchers::ContainsSubstring(format_significant(temperature) + " C"));
        // Inside its data nothing is logged.
        Core::get_remanence(material, lastRemanenceTemperature);
        CHECK(extrapolation_records(Logger::getInstance().drainCollected()).empty());
    }
    SECTION("missing data still throws with the flag on") {
        settings.set_allow_material_data_extrapolation(true);
        json emptyJson = materialJson;
        emptyJson["coerciveForce"] = json::array();
        CoreMaterial empty;
        from_json(emptyJson, empty);
        CHECK_THROWS_AS(Core::get_coercive_force(empty, temperature), InvalidInputException);
    }
    settings.reset();
}

namespace {

// The ErrorCode of what `call` throws as an OpenMagneticsException; none if it does not throw.
template <typename Call>
std::optional<ErrorCode> thrown_code(Call call) {
    try {
        call();
    }
    catch (const OpenMagneticsException& e) {
        return e.code();
    }
    return std::nullopt;
}

std::string thrown_message(const std::function<void()>& call) {
    try {
        call();
    }
    catch (const std::exception& e) {
        return e.what();
    }
    return "";
}

// 3C95 with its saturation, remanence, coercive force and resistivity replaced by the given tables, renamed so
// no cache keyed by name sees the catalogue's 3C95.
CoreMaterial material_with(const std::string& name, const std::function<void(json&)>& edit) {
    json materialJson;
    to_json(materialJson, Core::resolve_material("3C95"));
    materialJson["name"] = name;
    edit(materialJson);
    CoreMaterial material;
    from_json(materialJson, material);
    return material;
}

}  // namespace

// ABT #2235 (Alf, 2026-11-03): the saturation flux density, remanence, coercive force and saturation field strength
// below their first data point, a core material's resistivity and mu_i(f) outside their tabulated range are
// extrapolated exactly as they always were, whatever Settings::allowMaterialDataExtrapolation says: they never throw
// for being outside their data, and every such ask logs a WARNING from kMaterialDataExtrapolationModule. Inside the
// data nothing is logged. Where an extrapolated saturation flux density reaches zero or below it is unphysical and
// refused whatever the flag.
TEST_CASE("Saturation flux density outside its temperatures: both flags extrapolate and warn, unphysical always throws", "[material-extrapolation][saturation]") {
    settings.reset();
    clear_databases();
    CollectorGuard collector;
    auto material = material_with("3C95 with a two-point saturation table", [](json& j) {
        j["saturation"] = json::array({{{"magneticFluxDensity", 0.5}, {"magneticField", 1200}, {"temperature", 25}},
                                       {{"magneticFluxDensity", 0.4}, {"magneticField", 1200}, {"temperature", 100}}});
    });
    for (bool allow : {true, false}) {
        INFO("flag " << allow);
        settings.set_allow_material_data_extrapolation(allow);
        Logger::getInstance().drainCollected();
        CHECK(Core::get_magnetic_flux_density_saturation(material, 60, false) == Catch::Approx(0.5 - 35 * 0.1 / 75));
        CHECK(extrapolation_records(Logger::getInstance().drainCollected()).empty());
        double above = 0, below = 0;
        CHECK_NOTHROW(above = Core::get_magnetic_flux_density_saturation(material, 150, false));
        CHECK_NOTHROW(below = Core::get_magnetic_flux_density_saturation(material, -40, false));
        CHECK(above == Catch::Approx(0.5 - 125 * 0.1 / 75));
        CHECK(below == Catch::Approx(0.5 + 65 * 0.1 / 75));
        auto records = extrapolation_records(Logger::getInstance().drainCollected());
        REQUIRE(records.size() == 1);  // one material and quantity: the two asks merge into one record
        CHECK(records[0].level == LogLevel::WARNING);
        CHECK_THAT(records[0].message, Catch::Matchers::ContainsSubstring("3C95 with a two-point saturation table") &&
                                       Catch::Matchers::ContainsSubstring("saturation flux density") &&
                                       Catch::Matchers::ContainsSubstring("from 25 C to 100 C"));
        // 0.5 T - 475 C * 0.1 T / 75 C = -0.133 T at 500 C: unphysical, refused whatever the flag.
        CHECK(thrown_code([&] { Core::get_magnetic_flux_density_saturation(material, 500, false); }) == ErrorCode::MATERIAL_INVALID_PROPERTY);
        CHECK_THAT(thrown_message([&] { Core::get_magnetic_flux_density_saturation(material, 500, false); }),
                   Catch::Matchers::ContainsSubstring("is unphysical"));
    }
    settings.reset();
}

TEST_CASE("Saturation flux density from a one-point table: both flags hold the point and warn", "[material-extrapolation][saturation]") {
    settings.reset();
    clear_databases();
    CollectorGuard collector;
    auto material = material_with("3C95 with a one-point saturation table", [](json& j) {
        j["saturation"] = json::array({{{"magneticFluxDensity", 0.45}, {"magneticField", 1200}, {"temperature", 100}}});
    });
    for (bool allow : {true, false}) {
        INFO("flag " << allow);
        settings.set_allow_material_data_extrapolation(allow);
        Logger::getInstance().drainCollected();
        CHECK(Core::get_magnetic_flux_density_saturation(material, 100, false) == 0.45);
        CHECK(extrapolation_records(Logger::getInstance().drainCollected()).empty());
        double value = 0;
        CHECK_NOTHROW(value = Core::get_magnetic_flux_density_saturation(material, 20, false));
        CHECK(value == 0.45);
        auto records = extrapolation_records(Logger::getInstance().drainCollected());
        REQUIRE(records.size() == 1);
        CHECK_THAT(records[0].message, Catch::Matchers::ContainsSubstring("at 100 C only"));
        CHECK_THAT(records[0].detail, Catch::Matchers::ContainsSubstring("20 C"));
    }
    settings.reset();
}

// ABT #2235: below their first data point remanence, coercive force and the saturation field strength are their
// tables extrapolated (above it they run to zero at the Curie temperature, ABT #1736), whatever the flag, with a
// WARNING.
TEST_CASE("Ferromagnetic properties below their data: both flags extrapolate and warn", "[material-extrapolation][saturation]") {
    settings.reset();
    clear_databases();
    CollectorGuard collector;
    auto material = material_with("3C95 with two-point ferromagnetic tables", [](json& j) {
        j["saturation"] = json::array({{{"magneticFluxDensity", 0.5}, {"magneticField", 1200}, {"temperature", 25}},
                                       {{"magneticFluxDensity", 0.4}, {"magneticField", 1000}, {"temperature", 100}}});
        j["remanence"] = json::array({{{"magneticFluxDensity", 0.1}, {"magneticField", 0}, {"temperature", 25}},
                                      {{"magneticFluxDensity", 0.06}, {"magneticField", 0}, {"temperature", 100}}});
        j["coerciveForce"] = json::array({{{"magneticFluxDensity", 0}, {"magneticField", 15}, {"temperature", 25}},
                                          {{"magneticFluxDensity", 0}, {"magneticField", 10}, {"temperature", 100}}});
    });
    double temperature = -40;  // 65 C below the first point
    for (bool allow : {true, false}) {
        INFO("flag " << allow);
        settings.set_allow_material_data_extrapolation(allow);
        Logger::getInstance().drainCollected();
        CHECK(Core::get_remanence(material, 25) == Catch::Approx(0.1));
        CHECK(extrapolation_records(Logger::getInstance().drainCollected()).empty());
        double remanence = 0, coerciveForce = 0, saturationField = 0;
        CHECK_NOTHROW(remanence = Core::get_remanence(material, temperature));
        CHECK_NOTHROW(coerciveForce = Core::get_coercive_force(material, temperature));
        CHECK_NOTHROW(saturationField = Core::get_magnetic_field_strength_saturation(material, temperature));
        CHECK(remanence == Catch::Approx(0.1 + 65 * 0.04 / 75));
        CHECK(coerciveForce == Catch::Approx(15 + 65 * 5.0 / 75));
        CHECK(saturationField == Catch::Approx(1200 + 65 * 200.0 / 75));
        auto records = extrapolation_records(Logger::getInstance().drainCollected());
        REQUIRE(records.size() == 3);
        CHECK_THAT(records[0].message, Catch::Matchers::ContainsSubstring("remanence"));
        CHECK_THAT(records[1].message, Catch::Matchers::ContainsSubstring("coercive force"));
        CHECK_THAT(records[2].message, Catch::Matchers::ContainsSubstring("saturation field strength"));
        CHECK_THAT(records[0].detail, Catch::Matchers::ContainsSubstring("-40 C") && Catch::Matchers::ContainsSubstring("below"));
    }
    settings.reset();
}

// ABT #2235: a core material's resistivity outside its tabulated temperatures (one point, or several) is extrapolated
// as before (a single point is held), whatever the flag, with a WARNING.
TEST_CASE("Core resistivity outside its temperatures: both flags extrapolate and warn", "[material-extrapolation][resistivity]") {
    settings.reset();
    clear_databases();
    CollectorGuard collector;
    auto single = material_with("3C95 with one resistivity point", [](json& j) {
        j["resistivity"] = json::array({{{"value", 5}, {"temperature", 25}}});
    });
    // Collinear points: the table's spline and its extrapolation are that line.
    auto table = material_with("3C95 with a resistivity table", [](json& j) {
        j["resistivity"] = json::array({{{"value", 5}, {"temperature", 25}}, {{"value", 4}, {"temperature", 62.5}}, {{"value", 3}, {"temperature", 100}}});
    });
    for (bool allow : {true, false}) {
        INFO("flag " << allow);
        settings.set_allow_material_data_extrapolation(allow);
        Logger::getInstance().drainCollected();
        CHECK(Core::get_resistivity(single, 25) == 5);
        CHECK(Core::get_resistivity(table, 80) == Catch::Approx(5 - 55 * 2.0 / 75));
        CHECK(extrapolation_records(Logger::getInstance().drainCollected()).empty());
        double singleHot = 0, singleCold = 0, tableHot = 0, tableCold = 0;
        CHECK_NOTHROW(singleHot = Core::get_resistivity(single, 150));
        CHECK_NOTHROW(singleCold = Core::get_resistivity(single, -40));
        CHECK_NOTHROW(tableHot = Core::get_resistivity(table, 150));
        CHECK_NOTHROW(tableCold = Core::get_resistivity(table, -10));
        CHECK(singleHot == 5);
        CHECK(singleCold == 5);
        CHECK(tableHot == Catch::Approx(5 - 125 * 2.0 / 75));
        CHECK(tableCold == Catch::Approx(5 + 35 * 2.0 / 75));
        auto records = extrapolation_records(Logger::getInstance().drainCollected());
        REQUIRE(records.size() == 2);
        CHECK_THAT(records[0].message, Catch::Matchers::ContainsSubstring("3C95 with one resistivity point") &&
                                       Catch::Matchers::ContainsSubstring("resistivity") && Catch::Matchers::ContainsSubstring("at 25 C only"));
        CHECK_THAT(records[1].message, Catch::Matchers::ContainsSubstring("3C95 with a resistivity table") &&
                                       Catch::Matchers::ContainsSubstring("from 25 C to 100 C"));
    }
    settings.reset();
}

// ABT #2235: mu_i(f) outside its tabulated frequencies is the spline extrapolated (floored at 1) as before, whatever
// the flag, with a WARNING. PC95 tabulates mu_i from 2 kHz to 1.5 MHz.
TEST_CASE("Initial permeability outside its frequencies: both flags extrapolate the same and warn", "[material-extrapolation][initial-permeability]") {
    settings.reset();
    clear_databases();
    CollectorGuard collector;
    auto material = find_core_material_by_name("PC95");
    std::map<bool, std::vector<double>> values;
    for (bool allow : {true, false}) {
        INFO("flag " << allow);
        settings.set_allow_material_data_extrapolation(allow);
        Logger::getInstance().drainCollected();
        CHECK_NOTHROW(InitialPermeability::get_initial_permeability_frequency_dependent(material, 100e3));
        CHECK(extrapolation_records(Logger::getInstance().drainCollected()).empty());
        for (double frequency : {100.0, 500.0, 1.2e3, 2.9e6, 5e6, 20e6}) {
            double value = 0;
            CHECK_NOTHROW(value = InitialPermeability::get_initial_permeability_frequency_dependent(material, frequency));
            CHECK(value >= 1);
            values[allow].push_back(value);
        }
        auto records = extrapolation_records(Logger::getInstance().drainCollected());
        REQUIRE(records.size() == 1);
        CHECK(records[0].level == LogLevel::WARNING);
        CHECK_THAT(records[0].message, Catch::Matchers::ContainsSubstring("PC95") && Catch::Matchers::ContainsSubstring("initial permeability"));
    }
    CHECK(values[true] == values[false]);
    settings.reset();
}

// ABT #2235: the Magnetics temperature and frequency polynomials declare no span; outside the range they were
// fitted on they run below 1 (Edge 125 from 5 MHz). The checked getter refuses that whatever the flag, saying the
// value is unphysical. The raw model value (get_initial_permeability_model_value) stays the formula's value.
TEST_CASE("Magnetics initial permeability formula below 1: refused as unphysical whatever the flag", "[material-extrapolation][initial-permeability]") {
    settings.reset();
    clear_databases();
    auto material = find_core_material_by_name("Edge 125");
    for (bool allow : {true, false}) {
        INFO("flag " << allow);
        settings.set_allow_material_data_extrapolation(allow);
        CHECK_NOTHROW(InitialPermeability::get_initial_permeability(material, 25, std::nullopt, 100e3));
        CHECK_THROWS_AS(InitialPermeability::get_initial_permeability(material, 25, std::nullopt, 10e6), PermeabilityBelowVacuumException);
        CHECK_THAT(thrown_message([&] { InitialPermeability::get_initial_permeability(material, 25, std::nullopt, 10e6); }),
                   Catch::Matchers::ContainsSubstring("is unphysical"));
    }
    settings.reset();
}

// Alf, 2026-11-03: an adviser reads tables of materials the user never chose (a Vienna run logged 50,514 such
// WARNINGs, mostly Kool Mu saturation tabulated at 100 C only), and the web showed them all as warnings. While an
// adviser runs (a MaterialDataExtrapolationBarrier exists) a material read outside its tables is a DEBUG record;
// outside any adviser the same read is a WARNING, as ABT #2235 decided.
TEST_CASE("Material tables read outside their data under an adviser log DEBUG, outside it WARNING", "[material-extrapolation][adviser]") {
    settings.reset();
    clear_databases();
    CollectorGuard collector(LogLevel::DEBUG);
    auto material = material_with("3C95 with a saturation table at 100 C only", [](json& j) {
        j["saturation"] = json::array({{{"magneticFluxDensity", 0.45}, {"magneticField", 1200}, {"temperature", 100}}});
    });
    auto pc95 = find_core_material_by_name("PC95");  // mu_i tabulated from 2 kHz to 1.5 MHz
    auto levels = [&]() {
        CHECK(Core::get_magnetic_flux_density_saturation(material, 25, false) == 0.45);
        CHECK_NOTHROW(InitialPermeability::get_initial_permeability_frequency_dependent(pc95, 5e6));
        auto records = extrapolation_records(Logger::getInstance().drainCollected());
        REQUIRE(records.size() == 2);
        std::set<LogLevel> found;
        for (const auto& record : records) {
            found.insert(record.level);
        }
        return found;
    };
    Logger::getInstance().drainCollected();
    CHECK(levels() == std::set<LogLevel>{LogLevel::WARNING});
    {
        MaterialDataExtrapolationBarrier adviser;
        CHECK(levels() == std::set<LogLevel>{LogLevel::DEBUG});
    }
    CHECK(levels() == std::set<LogLevel>{LogLevel::WARNING});
    settings.reset();
}
