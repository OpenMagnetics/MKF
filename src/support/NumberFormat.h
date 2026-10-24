#pragma once
#include <string>

namespace OpenMagnetics {

// Numbers for the messages a person reads (findings, warnings): never std::to_string's six fixed
// decimals. A value keeps `significantFigures` significant figures, or all its integer digits when
// it has more (2474.39 -> "2474", 26.8475 -> "26.8", 0.25 -> "0.25"); trailing zeros are dropped.
// Values below 1e-4 in magnitude, and non-finite ones, are written in scientific notation.
std::string format_significant(double value, int significantFigures = 3);

// A quantity in an SI unit, scaled by the SI prefix (f, p, n, u, m, k, M, G, T) that puts its
// rounded mantissa in [1, 1000): (2.22e-5, "s") -> "22.2 us", (5.36e-4, "V*s") -> "536 uV*s",
// (0.027, "Ohm") -> "27 mOhm". Not for units that do not take prefixes (degrees C, ratios).
std::string format_quantity(double value, const std::string& unit, int significantFigures = 3);

} // namespace OpenMagnetics
