#include "support/NumberFormat.h"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace OpenMagnetics {

namespace {

// The decimal exponent of |value| (value != 0, finite): 2474 -> 3, 0.0536 -> -2.
int decimal_exponent(double value) {
    return static_cast<int>(std::floor(std::log10(std::fabs(value))));
}

// `value` rounded to `significantFigures` significant figures, with the digits it then has.
double round_significant(double value, int significantFigures) {
    int exponent = decimal_exponent(value);
    double scale = std::pow(10.0, exponent + 1 - significantFigures);
    return std::round(value / scale) * scale;
}

std::string drop_trailing_zeros(std::string text) {
    if (text.find('.') == std::string::npos) {
        return text;
    }
    while (!text.empty() && text.back() == '0') {
        text.pop_back();
    }
    if (!text.empty() && text.back() == '.') {
        text.pop_back();
    }
    return text;
}

} // namespace

std::string format_significant(double value, int significantFigures) {
    if (significantFigures < 1) {
        throw std::invalid_argument("format_significant: at least one significant figure is needed, not " +
                                    std::to_string(significantFigures));
    }
    std::ostringstream stream;
    if (!std::isfinite(value)) {
        stream << value;
        return stream.str();
    }
    if (value == 0) {
        return "0";
    }
    if (std::fabs(value) < 1e-4) {
        stream << std::setprecision(significantFigures) << value;
        return stream.str();
    }
    // Integer digits are never rounded away: 2474.39 C reads 2474 C, not 2470 C.
    int figures = std::max(significantFigures, decimal_exponent(value) + 1);
    double rounded = round_significant(value, figures);
    int decimals = std::max(0, figures - decimal_exponent(rounded) - 1);
    stream << std::fixed << std::setprecision(decimals) << rounded;
    std::string text = drop_trailing_zeros(stream.str());
    return text == "-0" ? "0" : text;
}

std::string format_quantity(double value, const std::string& unit, int significantFigures) {
    if (significantFigures < 1) {
        throw std::invalid_argument("format_quantity: at least one significant figure is needed, not " +
                                    std::to_string(significantFigures));
    }
    if (!std::isfinite(value) || value == 0) {
        return format_significant(value, significantFigures) + " " + unit;
    }
    static const char* prefixes[] = {"f", "p", "n", "u", "m", "", "k", "M", "G", "T"};
    constexpr int lowestPrefixPower = -15;
    constexpr int highestPrefixPower = 12;
    // The prefix is chosen on the rounded value, so 999.96 us reads 1 ms and not 1000 us.
    double rounded = round_significant(value, significantFigures);
    int exponent = decimal_exponent(rounded);
    int power = static_cast<int>(std::floor(exponent / 3.0)) * 3;
    if (power < lowestPrefixPower || power > highestPrefixPower) {
        std::ostringstream stream;
        stream << std::setprecision(significantFigures) << value << " " << unit;
        return stream.str();
    }
    std::string prefix = prefixes[(power - lowestPrefixPower) / 3];
    return format_significant(rounded / std::pow(10.0, power), significantFigures) + " " + prefix + unit;
}

} // namespace OpenMagnetics
