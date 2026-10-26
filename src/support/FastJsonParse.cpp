#include "support/FastJsonParse.h"
#include "simdjson.h"
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>

namespace OpenMagnetics {
namespace {

using json = nlohmann::json;

// simdjson reads the integer literal -0 as the integer 0, where nlohmann keeps a number_integer 0 (not the
// number_unsigned 0 of "0"). A text holding '-' directly followed by '0' and then not by '.', 'e', 'E' or a digit
// holds such a literal (or a string that looks like one); it goes to nlohmann.
bool holds_negative_zero_integer(const char* text, size_t length) {
    const char* end = text + length;
    for (const char* minus = static_cast<const char*>(std::memchr(text, '-', length)); minus;
         minus = static_cast<const char*>(std::memchr(minus + 1, '-', static_cast<size_t>(end - minus - 1)))) {
        if (minus + 1 >= end || minus[1] != '0') {
            continue;
        }
        const char* next = minus + 2;
        if (next >= end || !(*next == '.' || *next == 'e' || *next == 'E' || (*next >= '0' && *next <= '9'))) {
            return true;
        }
    }
    return false;
}

json convert(simdjson::dom::element element) {
    switch (element.type()) {
        case simdjson::dom::element_type::OBJECT: {
            json object = json::object();
            auto& members = object.get_ref<json::object_t&>();
            for (auto field : simdjson::dom::object(element)) {
                // operator[] then assignment: a repeated key keeps its last value, as nlohmann's parser does.
                members[std::string(field.key)] = convert(field.value);
            }
            return object;
        }
        case simdjson::dom::element_type::ARRAY: {
            json array = json::array();
            auto& items = array.get_ref<json::array_t&>();
            for (auto item : simdjson::dom::array(element)) {
                items.push_back(convert(item));
            }
            return array;
        }
        case simdjson::dom::element_type::INT64: {
            const int64_t value = int64_t(element);
            if (value >= 0) {
                return json(static_cast<json::number_unsigned_t>(value));
            }
            return json(static_cast<json::number_integer_t>(value));
        }
        case simdjson::dom::element_type::UINT64:
            return json(static_cast<json::number_unsigned_t>(uint64_t(element)));
        case simdjson::dom::element_type::DOUBLE:
            return json(static_cast<json::number_float_t>(double(element)));
        case simdjson::dom::element_type::STRING: {
            std::string_view value = std::string_view(element);
            return json(std::string(value.data(), value.size()));
        }
        case simdjson::dom::element_type::BOOL:
            return json(bool(element));
        case simdjson::dom::element_type::NULL_VALUE:
            return json(nullptr);
    }
    // Every simdjson element type is handled above; a new one must be added here, not guessed.
    throw std::logic_error("parse_json_fast: simdjson returned an element type this converter does not know");
}

} // namespace

nlohmann::json parse_json_fast(const char* text, size_t length) {
    if (!holds_negative_zero_integer(text, length)) {
        thread_local simdjson::dom::parser parser;
        simdjson::dom::element document;
        if (parser.parse(text, length).get(document) == simdjson::SUCCESS) {
            return convert(document);
        }
    }
    // Refused by simdjson, or a -0 literal: nlohmann's own result, or its exception.
    return nlohmann::json::parse(text, text + length);
}

} // namespace OpenMagnetics
