#pragma once
// Fast JSON parsing into the nlohmann::json tree MKF uses everywhere (Alf, 2026-10-26).
//
// simdjson parses the text; the result is converted into exactly the nlohmann::json value that
// nlohmann::json::parse would build from the same text: the same value types (a non-negative
// integer is number_unsigned, a negative one number_integer, a number with '.', 'e' or 'E'
// number_float), bit-identical doubles, sorted object keys with the last of duplicate keys
// kept, the same strings. Text simdjson refuses or cannot represent exactly (malformed JSON,
// integers beyond 64 bits, numbers out of double range, the integer literal -0) is parsed by
// nlohmann::json::parse itself, so a malformed document throws today's exception, text and all.
#include "json.hpp"
#include <cstddef>

namespace OpenMagnetics {

nlohmann::json parse_json_fast(const char* text, size_t length);

} // namespace OpenMagnetics
