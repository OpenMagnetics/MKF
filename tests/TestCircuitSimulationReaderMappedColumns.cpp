// ABT #2374 (web bug reports #196/#197, Alf's decisions 2026-11-13):
//  (a) only the columns the user maps are read as numbers. A NaN, an infinity or any other cell that is
//      not a number in a column nobody maps does not reject the file; one inside the period extracted
//      from a mapped column still throws, naming its line and column. Nothing is substituted.
//  (b) a file whose header has none of comma, semicolon, tab or pipe is split on runs of whitespace
//      (PSIM / SIMetrix .txt exports).
#include "processors/CircuitSimulatorInterface.h"
#include "support/Exceptions.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace OpenMagnetics;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;

namespace {

constexpr double kFrequency = 10e3;
constexpr double kTimeStep = 5e-7;
constexpr size_t kRows = 500;  // 2.5 periods of 10 kHz

// A 10 kHz triangle between 100 and 200 A.
double triangle(double t) {
    double phase = t * kFrequency - std::floor(t * kFrequency);
    return 100 + 100 * (phase < 0.5 ? 2 * phase : 2 * (1 - phase));
}

// One row per sample: time, current, voltage, then whatever extraCells(row) returns. The header is
// line 1, so row r is on line r + 2.
std::string build_csv(const std::string& header, const std::string& separator,
                      const std::function<std::string(size_t)>& currentCell,
                      const std::function<std::string(size_t)>& extraCells) {
    std::ostringstream csv;
    csv.precision(10);
    csv << header << "\n";
    for (size_t row = 0; row < kRows; ++row) {
        double t = row * kTimeStep;
        csv << t << separator << currentCell(row) << separator << 10;
        csv << extraCells(row) << "\n";
    }
    return csv.str();
}

std::string current_value(size_t row) {
    std::ostringstream cell;
    cell.precision(10);
    cell << triangle(row * kTimeStep);
    return cell.str();
}

const std::vector<std::map<std::string, std::string>> kMapping = {
    {{"time", "time"}, {"current", "I_L1"}, {"voltage", "V_L1"}}};

void check_imported_current(OperatingPoint operatingPoint) {
    auto current = operatingPoint.get_excitations_per_winding()[0].get_current().value();
    auto processed = current.get_processed().value();
    CHECK_THAT(processed.get_peak_to_peak().value(), WithinAbs(100, 2));
    CHECK_THAT(processed.get_offset(), WithinAbs(150, 2));
}

} // namespace

TEST_CASE("Reader: NaN, infinity and text in columns nobody maps do not reject the file",
          "[processor][circuit-simulation-reader][abt-2374][smoke-test]") {
    // P_ratio is NaN on every row (the X_nan_unused_col variant of the investigation); Z is an infinity
    // on the first row and text on another (Z_inf_unused).
    auto csv = build_csv("time,I_L1,V_L1,P_ratio,Z", ",", current_value, [](size_t row) {
        std::string z = row == 0 ? "inf" : (row == 300 ? "overflow" : "0");
        return std::string(",NaN,") + z;
    });
    CircuitSimulationReader reader(csv, true);
    REQUIRE(reader.extract_column_names().size() == 5);
    check_imported_current(reader.extract_operating_point(1, kFrequency, kMapping));

    // The automatic mapping the web offers first does not read those columns either.
    CircuitSimulationReader guessing(csv, true);
    auto guessed = guessing.extract_map_column_names(1, kFrequency);
    REQUIRE(guessed.size() == 1);
    CHECK(guessed[0].at("time") == "time");
    CHECK(guessed[0].at("current") == "I_L1");
}

TEST_CASE("Reader: a NaN in a mapped column inside the extracted period throws with its line and column",
          "[processor][circuit-simulation-reader][abt-2374][smoke-test]") {
    // Row 450 lies in the last period (rows 300..499), which is the one extracted.
    auto csv = build_csv("time,I_L1,V_L1", ",", [](size_t row) {
        return row == 450 ? std::string("nan") : current_value(row);
    }, [](size_t) { return std::string(); });
    CircuitSimulationReader reader(csv, true);
    try {
        reader.extract_operating_point(1, kFrequency, kMapping);
        FAIL("expected an exception");
    }
    catch (const InvalidInputException& e) {
        std::string message = e.what();
        CHECK_THAT(message, ContainsSubstring("line 452"));
        CHECK_THAT(message, ContainsSubstring("column 2"));
        CHECK_THAT(message, ContainsSubstring("\"I_L1\""));
        CHECK_THAT(message, ContainsSubstring("\"nan\""));
    }
}

TEST_CASE("Reader: a NaN in a mapped column before the extracted period does not reject the file",
          "[processor][circuit-simulation-reader][abt-2374][smoke-test]") {
    // The X2_nan_first_row variant: a simulator writes NaN on its first sample. Only the last period is
    // read, and the NaN is not in it.
    auto csv = build_csv("time,I_L1,V_L1", ",", [](size_t row) {
        return row == 0 ? std::string("NaN") : current_value(row);
    }, [](size_t) { return std::string(); });
    CircuitSimulationReader reader(csv, true);
    check_imported_current(reader.extract_operating_point(1, kFrequency, kMapping));
}

TEST_CASE("Reader: a time column with a cell that is not a number is refused, naming the cell",
          "[processor][circuit-simulation-reader][abt-2374][smoke-test]") {
    std::string csv = "time,I_L1\n0,1\nabc,3\n2e-6,2\n";
    try {
        CircuitSimulationReader reader(csv, true);
        FAIL("expected an exception");
    }
    catch (const InvalidInputException& e) {
        std::string message = e.what();
        CHECK_THAT(message, ContainsSubstring("no time column found"));
        CHECK_THAT(message, ContainsSubstring("line 3"));
        CHECK_THAT(message, ContainsSubstring("\"abc\""));
    }
}

TEST_CASE("guess_separator: runs of whitespace separate columns only when no other separator is present",
          "[processor][circuit-simulation-reader][abt-2374][smoke-test]") {
    CHECK(' ' == CircuitSimulationReader::guess_separator("Time I_L1 V_L1"));
    CHECK(' ' == CircuitSimulationReader::guess_separator(" Time         I_L1         V_L1"));
    CHECK(' ' == CircuitSimulationReader::guess_separator("time I(L1) V(p1,p2)"));
    // Names with spaces under a real separator keep that separator.
    CHECK(',' == CircuitSimulationReader::guess_separator("time,I(L1) A,V(out) V"));
    CHECK('\t' == CircuitSimulationReader::guess_separator("time\tI(L1) A\tV(out) V"));
    CHECK(';' == CircuitSimulationReader::guess_separator("Time / s;Current L1 / A"));
    REQUIRE_THROWS_AS(CircuitSimulationReader::guess_separator("   just_one_column   "), InvalidInputException);
}

TEST_CASE("Reader: a space-separated PSIM export imports with the right values",
          "[processor][circuit-simulation-reader][abt-2374][smoke-test]") {
    // E2_psim_multispace: a leading space on every line and several spaces between columns.
    std::ostringstream csv;
    csv.precision(7);
    csv << " Time         I_L1         V_L1\n";
    for (size_t row = 0; row < kRows; ++row) {
        double t = row * kTimeStep;
        csv << " " << std::scientific << t << "   " << triangle(t) << "   " << 24.0 << "\n";
    }
    CircuitSimulationReader reader(csv.str(), true);
    auto names = reader.extract_column_names();
    REQUIRE(names == std::vector<std::string>{"Time", "I_L1", "V_L1"});
    REQUIRE(reader.get_columns().at(1).data.size() == kRows);
    CHECK_THAT(reader.get_columns().at(1).data.at(0), WithinAbs(100, 1e-4));
    CHECK_THAT(reader.get_columns().at(2).data.at(kRows - 1), WithinAbs(24, 1e-9));
    check_imported_current(reader.extract_operating_point(
        1, kFrequency, std::vector<std::map<std::string, std::string>>{{{"time", "Time"}, {"current", "I_L1"}, {"voltage", "V_L1"}}}));
}
