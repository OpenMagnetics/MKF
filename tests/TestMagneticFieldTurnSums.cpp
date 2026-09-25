// The magnetic field's turn-field sums: summed once per mesh for every harmonic, and kept
// between calls (Settings::magnetic_field_turn_sums_cache_bytes). Keeping them must never change
// a result: every case compares against the same calculation with the cache off.
#include "physical_models/MagneticField.h"
#include "physical_models/WindingLosses.h"
#include "support/Settings.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace MAS;
using namespace OpenMagnetics;
using namespace OpenMagneticsTesting;

namespace {

// Restores the cache budget whatever a case does with it.
struct CacheBudgetGuard {
    size_t saved = settings.get_magnetic_field_turn_sums_cache_bytes();
    ~CacheBudgetGuard() {
        settings.set_magnetic_field_turn_sums_cache_bytes(saved);
        MagneticField::clear_turn_sums_cache();
    }
};

OperatingPoint triangular(double frequency, double peakToPeak, double dcCurrent, std::vector<double> turnsRatios = {}) {
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency, 100e-6, 25, WaveformLabel::TRIANGULAR, peakToPeak, 0.3, dcCurrent, turnsRatios);
    return inputs.get_operating_points()[0];
}

// A wound PQ 26/25 with a 0.5 mm ground gap.
OpenMagnetics::Magnetic part(std::vector<int64_t> numberTurns) {
    std::vector<int64_t> numberParallels(numberTurns.size(), 1);
    auto coil = get_quick_coil(numberTurns, numberParallels, "PQ 26/25");
    coil.delimit_and_compact();
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(get_quick_core("PQ 26/25", get_ground_gap(0.0005)));
    magnetic.set_coil(coil);
    return magnetic;
}

double winding_losses(OpenMagnetics::Magnetic magnetic, OperatingPoint operatingPoint) {
    return WindingLosses().calculate_losses(magnetic, operatingPoint, 25).get_winding_losses();
}

double uncached_winding_losses(OpenMagnetics::Magnetic magnetic, OperatingPoint operatingPoint) {
    size_t budget = settings.get_magnetic_field_turn_sums_cache_bytes();
    settings.set_magnetic_field_turn_sums_cache_bytes(0);
    double losses = winding_losses(magnetic, operatingPoint);
    settings.set_magnetic_field_turn_sums_cache_bytes(budget);
    return losses;
}

} // namespace

TEST_CASE("Turn-field sums are reused for another operating point of the same part", "[magnetic-field][turn-sums-cache]") {
    CacheBudgetGuard guard;
    MagneticField::clear_turn_sums_cache();
    auto magnetic = part({30});
    auto slow = triangular(100000, 2, 3);
    auto fast = triangular(400000, 0.5, 1);

    double first = winding_losses(magnetic, slow);
    double second = winding_losses(magnetic, fast);
    auto statistics = MagneticField::get_turn_sums_cache_statistics();

    // Summed once, then found: the second operating point reused the first one's sums.
    CHECK(statistics.misses == 1);
    CHECK(statistics.hits == 1);
    CHECK(statistics.entries == 1);
    CHECK_THAT(first, Catch::Matchers::WithinRel(uncached_winding_losses(magnetic, slow), 1e-12));
    CHECK_THAT(second, Catch::Matchers::WithinRel(uncached_winding_losses(magnetic, fast), 1e-12));
}

TEST_CASE("Turn-field sums of a two-winding part are reused, each winding with its own current", "[magnetic-field][turn-sums-cache]") {
    CacheBudgetGuard guard;
    MagneticField::clear_turn_sums_cache();
    auto magnetic = part({24, 12});
    auto light = triangular(100000, 1, 2, {2});
    auto heavy = triangular(200000, 3, 5, {2});

    double first = winding_losses(magnetic, light);
    double second = winding_losses(magnetic, heavy);

    CHECK(MagneticField::get_turn_sums_cache_statistics().hits == 1);
    CHECK_THAT(first, Catch::Matchers::WithinRel(uncached_winding_losses(magnetic, light), 1e-12));
    CHECK_THAT(second, Catch::Matchers::WithinRel(uncached_winding_losses(magnetic, heavy), 1e-12));
}

TEST_CASE("Another winding is another mesh: its turn-field sums are not reused", "[magnetic-field][turn-sums-cache]") {
    CacheBudgetGuard guard;
    MagneticField::clear_turn_sums_cache();
    auto operatingPoint = triangular(100000, 2, 3);
    auto thirty = part({30});
    auto twenty = part({20});

    winding_losses(thirty, operatingPoint);
    double losses = winding_losses(twenty, operatingPoint);

    auto statistics = MagneticField::get_turn_sums_cache_statistics();
    CHECK(statistics.hits == 0);
    CHECK(statistics.misses == 2);
    CHECK_THAT(losses, Catch::Matchers::WithinRel(uncached_winding_losses(twenty, operatingPoint), 1e-12));
}

TEST_CASE("The turn-field sums kept stay within the budget", "[magnetic-field][turn-sums-cache]") {
    CacheBudgetGuard guard;
    MagneticField::clear_turn_sums_cache();
    auto operatingPoint = triangular(100000, 2, 3);
    auto magnetic = part({30});

    // One part's sums, measured, then a budget that fits one part but not two.
    winding_losses(magnetic, operatingPoint);
    size_t onePart = MagneticField::get_turn_sums_cache_statistics().bytes;
    REQUIRE(onePart > 0);
    settings.set_magnetic_field_turn_sums_cache_bytes(onePart + onePart / 2);
    winding_losses(part({29}), operatingPoint);
    winding_losses(part({28}), operatingPoint);

    auto statistics = MagneticField::get_turn_sums_cache_statistics();
    CHECK(statistics.bytes <= onePart + onePart / 2);
    CHECK(statistics.entries == 1);

    // A budget of 0 keeps nothing.
    MagneticField::clear_turn_sums_cache();
    settings.set_magnetic_field_turn_sums_cache_bytes(0);
    winding_losses(magnetic, operatingPoint);
    CHECK(MagneticField::get_turn_sums_cache_statistics().entries == 0);
}
