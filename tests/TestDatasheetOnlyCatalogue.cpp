// =============================================================================
// TestDatasheetOnlyCatalogue.cpp
// =============================================================================
// A catalogue may hold parts known only by their datasheet: no core, no coil,
// just manufacturerInfo.datasheetInfo (MAS requires neither construction). The
// catalogue MagneticAdviser must rank them on what their datasheet states and
// must not fail, skip or invent a score for what it cannot judge:
//   - filters that need a construction do not apply (applies_to) and record no
//     score, so the part is ranked on the filters that do apply;
//   - the size filters read the published body size;
//   - MAGNETIZING_INDUCTANCE reads the measured L(I) at the operating bias and
//     temperature;
//   - DATASHEET_LIMITS reads the saturationCurrents table and the ΔT-qualified
//     rated currents, and scores utilisation (lower is better).
// =============================================================================

#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "advisers/MagneticAdviser.h"
#include "advisers/MagneticFilter.h"
#include "constructive_models/Magnetic.h"
#include "processors/Inputs.h"
#include "support/Exceptions.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "support/Cache.h"

#include "TestingUtils.h"

using namespace MAS;
using namespace OpenMagnetics;
using Catch::Matchers::WithinRel;

namespace {

DimensionWithTolerance nominal(double value) {
    DimensionWithTolerance dimension;
    dimension.set_nominal(value);
    return dimension;
}

struct DatasheetSpec {
    std::string reference;
    std::optional<double> inductance;
    // (current A, inductance H, temperature degC)
    std::vector<std::tuple<double, double, std::optional<double>>> curve;
    std::optional<double> ratedCurrent;
    std::vector<std::pair<double, double>> saturationCurrents;  // (percent drop, current)
    std::optional<std::tuple<double, double, double>> body;     // length, width, height (m)
};

// A datasheet-only single-winding inductor.
OpenMagnetics::Magnetic datasheet_part(const DatasheetSpec& spec) {
    MagneticDatasheetElectrical electrical;
    electrical.set_subtype(ElectricalSubtype::INDUCTOR);
    if (spec.inductance) {
        electrical.set_inductance(nominal(spec.inductance.value()));
    }
    if (!spec.curve.empty()) {
        std::vector<DatasheetInductancePoint> points;
        for (const auto& [current, inductance, temperature] : spec.curve) {
            DatasheetInductancePoint point;
            point.set_current(current);
            point.set_inductance(inductance);
            point.set_temperature(temperature);
            points.push_back(point);
        }
        electrical.set_inductance_points(points);
    }
    if (spec.ratedCurrent) {
        electrical.set_rated_currents(std::vector<double>{spec.ratedCurrent.value()});
    }
    if (!spec.saturationCurrents.empty()) {
        std::vector<DatasheetSaturationCurrent> table;
        for (const auto& [percent, current] : spec.saturationCurrents) {
            DatasheetSaturationCurrent entry;
            entry.set_percent_inductance_drop(percent);
            entry.set_current(current);
            table.push_back(entry);
        }
        electrical.set_saturation_currents(table);
    }
    DatasheetInfo datasheetInfo;
    datasheetInfo.set_electrical(std::vector<MagneticDatasheetElectrical>{electrical});
    if (spec.body) {
        Mechanical mechanical;
        mechanical.set_length(nominal(std::get<0>(spec.body.value())));
        mechanical.set_width(nominal(std::get<1>(spec.body.value())));
        mechanical.set_height(nominal(std::get<2>(spec.body.value())));
        datasheetInfo.set_mechanical(mechanical);
    }
    MagneticManufacturerInfo manufacturerInfo;
    manufacturerInfo.set_name("Test Maker");
    manufacturerInfo.set_reference(spec.reference);
    manufacturerInfo.set_datasheet_info(datasheetInfo);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_manufacturer_info(manufacturerInfo);
    return magnetic;
}

// 10 uH at zero bias, rolling off to 6 uH at 4 A at 20 C, and to 5 uH at 4 A at 100 C.
DatasheetSpec ten_microhenry(std::string reference) {
    DatasheetSpec spec;
    spec.reference = std::move(reference);
    spec.inductance = 10e-6;
    spec.curve = {{0, 10e-6, 20.0}, {2, 9e-6, 20.0}, {4, 6e-6, 20.0}, {0, 9e-6, 100.0}, {2, 8e-6, 100.0}, {4, 5e-6, 100.0}};
    spec.ratedCurrent = 3.0;
    spec.saturationCurrents = {{10, 2.5}, {30, 4.0}};
    spec.body = std::make_tuple(5e-3, 5e-3, 3e-3);
    return spec;
}

// 1 MHz triangular current around `dcCurrent` with 0.4 A ripple, asking for `inductance`.
OpenMagnetics::Inputs buck_inductor_inputs(double inductance, double dcCurrent, double temperature = 25) {
    return OpenMagnetics::Inputs::create_quick_operating_point_only_current(1e6, inductance, temperature, WaveformLabel::TRIANGULAR, 0.4, 0.5, dcCurrent);
}

}  // namespace

TEST_CASE("A datasheet-only part reads its inductance from the measured L(I)",
          "[datasheet-only][magnetic][smoke-test]") {
    auto magnetic = datasheet_part(ten_microhenry("P1"));

    SECTION("linear in current along a measured temperature") {
        REQUIRE_THAT(magnetic.calculate_datasheet_inductance(1.0, 20).value(), WithinRel(9.5e-6, 1e-9));
        REQUIRE_THAT(magnetic.calculate_datasheet_inductance(-3.0, 20).value(), WithinRel(7.5e-6, 1e-9));
    }
    SECTION("linear between the two bracketing temperatures") {
        REQUIRE_THAT(magnetic.calculate_datasheet_inductance(2.0, 60).value(), WithinRel(8.5e-6, 1e-9));
    }
    SECTION("the nearest measured temperature outside the measured range") {
        REQUIRE_THAT(magnetic.calculate_datasheet_inductance(2.0, -40).value(), WithinRel(9e-6, 1e-9));
        REQUIRE_THAT(magnetic.calculate_datasheet_inductance(2.0, 150).value(), WithinRel(8e-6, 1e-9));
    }
    SECTION("nothing beyond the last measured current") {
        REQUIRE_FALSE(magnetic.calculate_datasheet_inductance(4.5, 20).has_value());
    }
}

TEST_CASE("A datasheet-only part without L(I) reads its stated inductance, and one with neither throws",
          "[datasheet-only][magnetic][smoke-test]") {
    DatasheetSpec spec;
    spec.reference = "P2";
    spec.inductance = 4.7e-6;
    REQUIRE_THAT(datasheet_part(spec).calculate_datasheet_inductance(1.0, 25).value(), WithinRel(4.7e-6, 1e-9));

    spec.inductance.reset();
    spec.ratedCurrent = 1.0;
    REQUIRE_THROWS_AS(datasheet_part(spec).calculate_datasheet_inductance(1.0, 25), InvalidInputException);
}

TEST_CASE("L(I) points with and without a temperature are refused, not read one way or the other",
          "[datasheet-only][magnetic][smoke-test]") {
    DatasheetSpec spec;
    spec.reference = "P3";
    spec.curve = {{0, 10e-6, 20.0}, {2, 9e-6, std::nullopt}};
    REQUIRE_THROWS_AS(datasheet_part(spec).calculate_datasheet_inductance(1.0, 25), InvalidInputException);
}

TEST_CASE("A datasheet-only part's dimensions are its published body size",
          "[datasheet-only][magnetic][smoke-test]") {
    auto magnetic = datasheet_part(ten_microhenry("P4"));
    REQUIRE(magnetic.has_datasheet_dimensions());
    auto dimensions = magnetic.get_maximum_dimensions();
    REQUIRE(dimensions == std::vector<double>{5e-3, 3e-3, 5e-3});

    DatasheetSpec bare;
    bare.reference = "P5";
    bare.inductance = 1e-6;
    auto sizeless = datasheet_part(bare);
    REQUIRE_FALSE(sizeless.has_datasheet_dimensions());
    REQUIRE_THROWS_AS(sizeless.get_maximum_dimensions(), InvalidInputException);
}

TEST_CASE("Filters that need a construction do not apply to a datasheet-only part; datasheet filters do",
          "[datasheet-only][magnetic-filter][smoke-test]") {
    auto magnetic = datasheet_part(ten_microhenry("P6"));
    for (auto filter : {MagneticFilters::LOSSES, MagneticFilters::SATURATION, MagneticFilters::TEMPERATURE_RISE, MagneticFilters::DC_CURRENT_DENSITY}) {
        INFO(std::string(magic_enum::enum_name(filter)));
        REQUIRE_FALSE(MagneticFilter::factory(filter)->applies_to(&magnetic));
    }
    for (auto filter : {MagneticFilters::DATASHEET_LIMITS, MagneticFilters::MAGNETIZING_INDUCTANCE, MagneticFilters::VOLUME, MagneticFilters::AREA, MagneticFilters::HEIGHT, MagneticFilters::MAXIMUM_DIMENSIONS}) {
        INFO(std::string(magic_enum::enum_name(filter)));
        REQUIRE(MagneticFilter::factory(filter)->applies_to(&magnetic));
    }
}

TEST_CASE("MAGNETIZING_INDUCTANCE judges a datasheet-only part at its operating bias",
          "[datasheet-only][magnetic-filter][smoke-test]") {
    auto magnetic = datasheet_part(ten_microhenry("P7"));
    auto filter = MagneticFilter::factory(MagneticFilters::MAGNETIZING_INDUCTANCE);

    SECTION("enough inductance left at 1 A") {
        auto inputs = buck_inductor_inputs(9e-6, 1.0, 20);
        auto [valid, score] = filter->evaluate_magnetic(&magnetic, &inputs);
        REQUIRE(valid);
        REQUIRE_THAT(score, WithinRel(0.5e-6, 1e-6));
    }
    SECTION("too little left at 3.5 A, where the curve has rolled off") {
        auto inputs = buck_inductor_inputs(9e-6, 3.5, 20);
        auto [valid, score] = filter->evaluate_magnetic(&magnetic, &inputs);
        REQUIRE_FALSE(valid);
    }
    SECTION("a bias beyond the measured curve fails") {
        auto inputs = buck_inductor_inputs(1e-6, 5.0, 20);
        auto [valid, score] = filter->evaluate_magnetic(&magnetic, &inputs);
        REQUIRE_FALSE(valid);
    }
}

TEST_CASE("DATASHEET_LIMITS gates the peak against the tightest saturation criterion the datasheet states",
          "[datasheet-only][datasheet-limits][smoke-test]") {
    // Isat is 2.5 A at 10 % and 4.0 A at 30 %: a 2.8 A peak passes the 30 % figure but not the 10 % one.
    auto magnetic = datasheet_part(ten_microhenry("P8"));
    auto filter = MagneticFilter::factory(MagneticFilters::DATASHEET_LIMITS);
    auto inputs = buck_inductor_inputs(10e-6, 2.6);  // peak 2.8 A, RMS below the 3 A rating
    auto [valid, score] = filter->evaluate_magnetic(&magnetic, &inputs);
    REQUIRE_FALSE(valid);
    REQUIRE(score > 1.0);
}

TEST_CASE("DATASHEET_LIMITS uses the smallest ΔT-qualified rated current when no plain rating is given",
          "[datasheet-only][datasheet-limits][smoke-test]") {
    MagneticDatasheetElectrical electrical;
    electrical.set_subtype(ElectricalSubtype::INDUCTOR);
    DatasheetRatedCurrent at20;
    at20.set_current(2.0);
    at20.set_temperature_rise(20);
    DatasheetRatedCurrent at40;
    at40.set_current(3.0);
    at40.set_temperature_rise(40);
    electrical.set_rated_current_points(std::vector<DatasheetRatedCurrent>{at20, at40});
    auto magnetic = datasheet_part(DatasheetSpec{"P9"});
    auto manufacturerInfo = magnetic.get_manufacturer_info().value();
    auto datasheetInfo = manufacturerInfo.get_datasheet_info().value();
    datasheetInfo.set_electrical(std::vector<MagneticDatasheetElectrical>{electrical});
    manufacturerInfo.set_datasheet_info(datasheetInfo);
    magnetic.set_manufacturer_info(manufacturerInfo);

    auto filter = MagneticFilter::factory(MagneticFilters::DATASHEET_LIMITS);
    auto inputs = buck_inductor_inputs(10e-6, 2.5);
    auto [valid, score] = filter->evaluate_magnetic(&magnetic, &inputs);
    REQUIRE_FALSE(valid);
    REQUIRE_THAT(score, WithinRel(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_rms().value() / 2.0, 1e-9));
}

TEST_CASE("DATASHEET_LIMITS ranks the part with more margin first under the usual invert=true",
          "[datasheet-only][datasheet-limits][smoke-test]") {
    // Utilisation is lower-is-better, like every other filter. Before 2026-09 the filter returned
    // headroom (higher-is-better) and invert=true ranked the part closest to its ratings first.
    auto tight = ten_microhenry("TIGHT");
    tight.ratedCurrent = 1.2;
    tight.saturationCurrents = {{30, 10.0}};
    auto roomy = ten_microhenry("ROOMY");
    roomy.ratedCurrent = 5.0;
    roomy.saturationCurrents = {{30, 10.0}};
    std::vector<OpenMagnetics::Magnetic> catalogue{datasheet_part(tight), datasheet_part(roomy)};

    std::vector<MagneticFilterOperation> flow{MagneticFilterOperation(MagneticFilters::DATASHEET_LIMITS, true, false, false, 1.0)};
    MagneticAdviser adviser;
    auto results = adviser.get_advised_magnetic(buck_inductor_inputs(10e-6, 1.0), catalogue, flow, 2, false);
    REQUIRE(results.size() == 2);
    REQUIRE(results[0].first.get_magnetic().get_reference() == "ROOMY");
}

TEST_CASE("The catalogue adviser ranks datasheet-only parts next to constructed ones",
          "[datasheet-only][adviser]") {
    settings.reset();
    OpenMagneticsTesting::QuickMagneticConfig cfg;
    cfg.numberTurns = {10};
    cfg.numberParallels = {1};
    cfg.coreShapeName = "E 35";
    cfg.coreMaterialName = "3C97";
    cfg.wireNames = {"Round 1.00 - Grade 1"};
    cfg.numberStacks = 1;
    auto constructed = OpenMagneticsTesting::create_quick_test_magnetic(cfg);
    constructed.get_mutable_coil().wind();
    MagneticManufacturerInfo constructedInfo;
    constructedInfo.set_name("Test Maker");
    constructedInfo.set_reference("BUILT");
    constructed.set_manufacturer_info(constructedInfo);

    std::vector<OpenMagnetics::Magnetic> catalogue{constructed, datasheet_part(ten_microhenry("SHEET-A")), datasheet_part(ten_microhenry("SHEET-B"))};
    std::vector<MagneticFilterOperation> flow{
        MagneticFilterOperation(MagneticFilters::DATASHEET_LIMITS, true, false, true, 1.0),
        MagneticFilterOperation(MagneticFilters::MAGNETIZING_INDUCTANCE, true, false, false, 1.0),
        MagneticFilterOperation(MagneticFilters::VOLUME, true, false, false, 1.0),
        MagneticFilterOperation(MagneticFilters::LOSSES, true, false, false, 1.0),
    };
    MagneticAdviser adviser;
    std::vector<std::pair<OpenMagnetics::Mas, double>> results;
    REQUIRE_NOTHROW(results = adviser.get_advised_magnetic(buck_inductor_inputs(9e-6, 1.0), catalogue, flow, 10, false));

    std::set<std::string> references;
    for (auto& [mas, score] : results) {
        references.insert(mas.get_magnetic().get_reference());
        CHECK(std::isfinite(score));
    }
    // Both datasheet parts survive to the results -- including the final simulation, which has
    // nothing to simulate for them and must not drop them.
    REQUIRE(references.count("SHEET-A") == 1);
    REQUIRE(references.count("SHEET-B") == 1);

    auto scorings = adviser.get_scorings();
    REQUIRE(scorings["SHEET-A"].count(MagneticFilters::MAGNETIZING_INDUCTANCE) == 1);
    REQUIRE(scorings["SHEET-A"].count(MagneticFilters::VOLUME) == 1);
    // LOSSES cannot judge a part with no construction: it records nothing, never a stand-in.
    REQUIRE(scorings["SHEET-A"].count(MagneticFilters::LOSSES) == 0);
    REQUIRE(scorings["BUILT"].count(MagneticFilters::LOSSES) == 1);
    settings.reset();
}

TEST_CASE("A datasheet-only part whose winding count cannot be told is refused loudly",
          "[datasheet-only][adviser][smoke-test]") {
    auto magnetic = datasheet_part(ten_microhenry("P10"));
    auto manufacturerInfo = magnetic.get_manufacturer_info().value();
    auto datasheetInfo = manufacturerInfo.get_datasheet_info().value();
    auto electrical = datasheetInfo.get_electrical().value();
    electrical[0].set_subtype(ElectricalSubtype::COUPLED_INDUCTOR);
    datasheetInfo.set_electrical(electrical);
    manufacturerInfo.set_datasheet_info(datasheetInfo);
    magnetic.set_manufacturer_info(manufacturerInfo);

    MagneticAdviser adviser;
    std::vector<MagneticFilterOperation> flow{MagneticFilterOperation(MagneticFilters::DATASHEET_LIMITS, true, false, false, 1.0)};
    REQUIRE_THROWS_AS(adviser.get_advised_magnetic(buck_inductor_inputs(10e-6, 1.0), std::vector<OpenMagnetics::Magnetic>{magnetic}, flow, 1, false), InvalidInputException);
}

TEST_CASE("Autocompleting a datasheet-only part leaves it as it is",
          "[datasheet-only][magnetic][smoke-test]") {
    auto magnetic = datasheet_part(ten_microhenry("P11"));
    OpenMagnetics::Magnetic completed;
    REQUIRE_NOTHROW(completed = magnetic_autocomplete(magnetic));
    REQUIRE_FALSE(completed.has_core());
    REQUIRE(completed.get_reference() == "P11");
}

TEST_CASE("A cache subset holds exactly the requested references, and refuses unknown ones",
          "[datasheet-only][cache][smoke-test]") {
    Cache<OpenMagnetics::Magnetic> cache;
    for (auto reference : {"A", "B", "C"}) {
        cache.load(reference, datasheet_part(ten_microhenry(reference)));
    }
    auto subset = cache.subset({"A", "C"});
    REQUIRE(subset.size() == 2);
    REQUIRE(subset.count("A") == 1);
    REQUIRE(subset.count("B") == 0);
    REQUIRE_THROWS_WITH(cache.subset({"A", "Z"}), Catch::Matchers::ContainsSubstring("1 of 2 requested references are not in the cache (Z)"));
}

TEST_CASE("The catalogue adviser weights each filter once: the score is the weighted mean",
          "[datasheet-only][adviser][scoring][smoke-test]") {
    // Three parts, two filters: VOLUME (weight 2) and DATASHEET_LIMITS (weight 1). SMALL is the
    // smallest and the most loaded; BIG the largest with the most margin; MID sits between.
    // Normalised per filter (best 1, worst 0), SMALL scores VOLUME 1 and DATASHEET_LIMITS 0, so
    // its weighted mean is (2*1 + 1*0) / 3. Weighting twice gave (2*2 + 0) / 3 = 1.33.
    auto small = ten_microhenry("SMALL");
    small.ratedCurrent = 1.2;
    small.body = std::make_tuple(3e-3, 3e-3, 2e-3);
    auto mid = ten_microhenry("MID");
    mid.ratedCurrent = 2.0;
    mid.body = std::make_tuple(4e-3, 4e-3, 3e-3);
    auto big = ten_microhenry("BIG");
    big.ratedCurrent = 6.0;
    big.body = std::make_tuple(6e-3, 6e-3, 4e-3);
    for (auto* spec : {&small, &mid, &big}) {
        spec->saturationCurrents = {{30, 20.0}};
    }
    std::vector<OpenMagnetics::Magnetic> catalogue{datasheet_part(small), datasheet_part(mid), datasheet_part(big)};
    std::vector<MagneticFilterOperation> flow{
        MagneticFilterOperation(MagneticFilters::VOLUME, true, false, false, 2.0),
        MagneticFilterOperation(MagneticFilters::DATASHEET_LIMITS, true, false, false, 1.0),
    };
    MagneticAdviser adviser;
    auto results = adviser.get_advised_magnetic(buck_inductor_inputs(10e-6, 1.0), catalogue, flow, 3, false);
    REQUIRE(results.size() == 3);
    std::map<std::string, double> score;
    for (auto& [mas, value] : results) {
        score[mas.get_magnetic().get_reference()] = value;
        REQUIRE(value >= 0.0);
        REQUIRE(value <= 1.0);
    }
    REQUIRE_THAT(score["SMALL"], WithinRel(2.0 / 3.0, 1e-9));
    REQUIRE_THAT(score["BIG"], WithinRel(1.0 / 3.0, 1e-9));
}

namespace {

// A winding current over one period at 500 kHz: `samples` evenly spread, the last repeating the first.
OperatingPointExcitation winding_current(const std::vector<double>& samples) {
    const double frequency = 500000;
    std::vector<double> time;
    for (size_t i = 0; i < samples.size(); ++i) {
        time.push_back(i / (frequency * (samples.size() - 1)));
    }
    Waveform waveform;
    waveform.set_data(samples);
    waveform.set_time(time);
    SignalDescriptor current;
    current.set_waveform(waveform);
    current.set_processed(OpenMagnetics::Inputs::calculate_processed_data(waveform, frequency));
    OperatingPointExcitation excitation;
    excitation.set_frequency(frequency);
    excitation.set_current(current);
    return excitation;
}

OpenMagnetics::Inputs coupled_inputs(double inductance, const std::vector<std::vector<double>>& windings,
                                     std::vector<IsolationSide> sides = {IsolationSide::PRIMARY, IsolationSide::PRIMARY},
                                     std::vector<double> turnsRatios = {}) {
    OpenMagnetics::Inputs inputs;
    DimensionWithTolerance required;
    required.set_minimum(inductance);
    inputs.get_mutable_design_requirements().set_magnetizing_inductance(required);
    std::vector<DimensionWithTolerance> ratios;
    for (double ratio : turnsRatios) {
        ratios.push_back(nominal(ratio));
    }
    inputs.get_mutable_design_requirements().set_turns_ratios(ratios);
    inputs.get_mutable_design_requirements().set_isolation_sides(sides);
    OperatingPoint operatingPoint;
    OperatingConditions conditions;
    conditions.set_ambient_temperature(20);
    operatingPoint.set_conditions(conditions);
    for (const auto& samples : windings) {
        operatingPoint.get_mutable_excitations_per_winding().push_back(winding_current(samples));
    }
    inputs.get_mutable_operating_points().push_back(operatingPoint);
    return inputs;
}

// The coupled inductor as Heimdall stores it: a 'single winding' inductor entry (with L(I)) and a
// coupledInductor entry restating the per-winding values, two windings.
OpenMagnetics::Magnetic coupled_part(const std::string& reference, bool coupled) {
    auto magnetic = datasheet_part(ten_microhenry(reference));
    auto manufacturerInfo = magnetic.get_manufacturer_info().value();
    auto datasheetInfo = manufacturerInfo.get_datasheet_info().value();
    auto electrical = datasheetInfo.get_electrical().value();
    if (coupled) {
        auto pair = electrical[0];
        pair.set_subtype(ElectricalSubtype::COUPLED_INDUCTOR);
        pair.set_inductance_points(std::nullopt);
        electrical.push_back(pair);
    }
    datasheetInfo.set_electrical(electrical);
    Part part;
    part.set_number_of_windings(2);
    datasheetInfo.set_part(part);
    manufacturerInfo.set_datasheet_info(datasheetInfo);
    magnetic.set_manufacturer_info(manufacturerInfo);
    return magnetic;
}

const std::vector<double> kRising = {0.9, 1.2, 1.5, 1.2, 0.9};   // DC 1.2 A, peak 1.5 A
const std::vector<double> kFalling = {1.5, 1.2, 0.9, 1.2, 1.5};  // the same, half a period later

double ampere_turn_peak(const OpenMagnetics::Inputs& inputs, std::vector<IsolationSide> sides, std::vector<double> turnsRatios = {}) {
    return OpenMagnetics::Inputs::calculate_ampere_turn_current(inputs.get_operating_points()[0], turnsRatios, sides).get_processed()->get_peak().value();
}

}  // namespace

TEST_CASE("The ampere-turn current sums the windings in time, with phase, side and turns",
          "[datasheet-only][coupled][smoke-test]") {
    const std::vector<IsolationSide> primaries{IsolationSide::PRIMARY, IsolationSide::PRIMARY};
    SECTION("in phase, two 1.5 A peaks make a 3 A peak") {
        REQUIRE_THAT(ampere_turn_peak(coupled_inputs(5e-6, {kRising, kRising}), primaries), WithinRel(3.0, 1e-3));
    }
    SECTION("in anti-phase the equal ripples cancel: only the two 1.2 A DC components remain, 2.4 A") {
        // kRising and kFalling are 1.2 A DC with +/-0.3 A of ripple half a period apart, so the
        // ripple sums to zero at every sample and the sum is a flat 2.4 A -- the DC part alone.
        REQUIRE_THAT(ampere_turn_peak(coupled_inputs(5e-6, {kRising, kFalling}), primaries), WithinRel(2.4, 1e-3));
    }
    SECTION("pure AC of equal amplitude in anti-phase cancels completely") {
        const std::vector<double> ac = {-0.3, 0.0, 0.3, 0.0, -0.3};
        const std::vector<double> antiAc = {0.3, 0.0, -0.3, 0.0, 0.3};
        REQUIRE(ampere_turn_peak(coupled_inputs(5e-6, {ac, antiAc}), primaries) < 1e-9);
    }
    SECTION("a secondary-side winding counts negative") {
        std::vector<IsolationSide> sides{IsolationSide::PRIMARY, IsolationSide::SECONDARY};
        REQUIRE(ampere_turn_peak(coupled_inputs(5e-6, {kRising, kRising}, sides), sides) < 1e-9);
    }
    SECTION("a winding of half the turns (Np/Ns = 2) counts half") {
        REQUIRE_THAT(ampere_turn_peak(coupled_inputs(5e-6, {kRising, kRising}, primaries, {2.0}), primaries, {2.0}), WithinRel(2.25, 1e-3));
    }
    SECTION("processed values alone carry no phase: a winding without a waveform throws") {
        auto inputs = coupled_inputs(5e-6, {kRising, kRising});
        auto current = inputs.get_mutable_operating_points()[0].get_mutable_excitations_per_winding()[1].get_current().value();
        current.set_waveform(std::nullopt);
        inputs.get_mutable_operating_points()[0].get_mutable_excitations_per_winding()[1].set_current(current);
        REQUIRE_THROWS_AS(OpenMagnetics::Inputs::calculate_ampere_turn_current(inputs.get_operating_points()[0], {}, primaries), InvalidInputException);
    }
}

TEST_CASE("A datasheet coupled inductor is gated by the peak of its ampere-turn current",
          "[datasheet-only][coupled][smoke-test]") {
    // Isat 2.5 A at 10 % per winding; each winding peaks at 1.5 A (under its 3 A rating).
    auto filter = MagneticFilter::factory(MagneticFilters::DATASHEET_LIMITS);
    auto coupled = coupled_part("COUPLED", true);
    REQUIRE(coupled.is_datasheet_coupled_inductor());

    SECTION("in phase: 3 A of ampere-turn current saturates it") {
        auto inputs = coupled_inputs(5e-6, {kRising, kRising});
        auto [valid, score] = filter->evaluate_magnetic(&coupled, &inputs);
        REQUIRE_FALSE(valid);
        REQUIRE_THAT(score, WithinRel(3.0 / 2.5, 1e-3));
    }
    SECTION("in anti-phase: 2.4 A does not -- summing the peaks (3 A) would have rejected it") {
        auto inputs = coupled_inputs(5e-6, {kRising, kFalling});
        auto [valid, score] = filter->evaluate_magnetic(&coupled, &inputs);
        REQUIRE(valid);
    }
    SECTION("a two-winding datasheet part NOT stated as coupled keeps the per-winding peak") {
        auto uncoupled = coupled_part("TWO", false);
        REQUIRE_FALSE(uncoupled.is_datasheet_coupled_inductor());
        auto inputs = coupled_inputs(5e-6, {kRising, kRising});
        auto [valid, score] = filter->evaluate_magnetic(&uncoupled, &inputs);
        REQUIRE(valid);
    }
    SECTION("without isolation sides the windings cannot be signed: it throws") {
        auto inputs = coupled_inputs(5e-6, {kRising, kRising});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::nullopt);
        REQUIRE_THROWS_AS(filter->evaluate_magnetic(&coupled, &inputs), InvalidInputException);
    }
}

TEST_CASE("A datasheet coupled inductor's inductance is read at its ampere-turn DC",
          "[datasheet-only][coupled][smoke-test]") {
    // L(I) at 20 C: 10 uH at 0 A, 9 uH at 2 A, 6 uH at 4 A. 1.5 A DC in each winding, in phase:
    // 3 A of ampere-turn DC, 7.5 uH -- below the 8 uH asked for (1.5 A alone would leave 9.25 uH).
    const std::vector<double> dc15 = {1.2, 1.5, 1.8, 1.5, 1.2};
    auto magnetic = coupled_part("COUPLED", true);
    auto inputs = coupled_inputs(8e-6, {dc15, dc15});
    auto filter = MagneticFilter::factory(MagneticFilters::MAGNETIZING_INDUCTANCE);
    auto [valid, score] = filter->evaluate_magnetic(&magnetic, &inputs);
    REQUIRE_FALSE(valid);
    REQUIRE_THAT(score, WithinRel(0.5e-6, 1e-3));
}

TEST_CASE("A two-winding operating point searches the coupled inductors of a catalogue",
          "[datasheet-only][coupled][adviser][smoke-test]") {
    std::vector<OpenMagnetics::Magnetic> catalogue{coupled_part("COUPLED", true), datasheet_part(ten_microhenry("SINGLE"))};
    std::vector<MagneticFilterOperation> flow{MagneticFilterOperation(MagneticFilters::DATASHEET_LIMITS, true, false, false, 1.0)};
    MagneticAdviser adviser;
    const std::vector<double> small = {0.3, 0.5, 0.7, 0.5, 0.3};
    auto results = adviser.get_advised_magnetic(coupled_inputs(5e-6, {small, small}), catalogue, flow, 5, false);
    REQUIRE(results.size() == 1);
    REQUIRE(results[0].first.get_magnetic().get_reference() == "COUPLED");
}
