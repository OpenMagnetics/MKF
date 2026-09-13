#pragma once
// ABT #1210: how a family bobbin processor hands its geometry to BobbinDataProcessor.
//
// The generated CoreBobbinProcessedDescription holds column_depth, column_shape, column_thickness
// and wall_thickness as PLAIN members with no initialiser, so a processor that forgot one returned
// an indeterminate value -- undefined behaviour, not an exception. The E and EFD processors (and
// ER/EL/P/U, which share the E one) never set column_depth, and every turn wound on those bobbins
// came out short by four times the column half-depth. A family processor therefore no longer builds
// the description itself: it fills a BobbinFamilyGeometry, whose members are all optionals, and
// assemble_bobbin_processed_description refuses (throws, naming the member) to build a description
// from one that leaves any of them unset.
//
// Deliberately kept out of Bobbin.h: only Bobbin.cpp and its tests need it.

#include "constructive_models/Bobbin.h"
#include "support/Exceptions.h"
#include "support/Utils.h"

#include <map>
#include <optional>
#include <string>
#include <vector>

namespace OpenMagnetics {

struct BobbinFamilyGeometry {
    std::optional<ColumnShape> columnShape;
    std::optional<double> columnThickness;
    std::optional<double> wallThickness;
    // Half-extents of the central column INCLUDING the bobbin's column wall, the convention
    // create_quick_bobbin uses (core half-dimension + column thickness). A round column's are both
    // its tube radius.
    std::optional<double> columnWidth;
    std::optional<double> columnDepth;
    std::optional<WindingWindowElement> windingWindow;
};

// A bobbin's labelled dimensions. Asking for a label the record does not declare throws, naming the
// label, the bobbin and the labels it does declare (flatten_dimensions alone answers 0 for a missing
// key through map::operator[]).
class BobbinLabelledDimensions {
  public:
    explicit BobbinLabelledDimensions(OpenMagnetics::Bobbin& bobbin) {
        if (!bobbin.get_functional_description()) {
            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
                "Bobbin has no functionalDescription, so it has no labelled dimensions to process.");
        }
        auto functionalDescription = bobbin.get_functional_description().value();
        _values = flatten_dimensions(functionalDescription.get_dimensions());
        _bobbinName = bobbin.get_name() ? bobbin.get_name().value() : std::string("<unnamed>");
        _familyName = to_string(functionalDescription.get_family());
    }

    double operator()(const std::string& label) const {
        auto it = _values.find(label);
        if (it == _values.end()) {
            std::string declared;
            for (const auto& [key, _] : _values) {
                if (!declared.empty()) {
                    declared += ", ";
                }
                declared += key;
            }
            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
                "Bobbin '" + _bobbinName + "' (family '" + _familyName + "') has no dimension '" + label +
                "', which its family processor needs; the record declares {" + declared + "}.");
        }
        return it->second;
    }

    const std::string& bobbin_name() const { return _bobbinName; }
    const std::string& family_name() const { return _familyName; }

  private:
    std::map<std::string, double> _values;
    std::string _bobbinName;
    std::string _familyName;
};

inline CoreBobbinProcessedDescription assemble_bobbin_processed_description(const BobbinFamilyGeometry& geometry,
                                                                            const std::string& familyName,
                                                                            const std::string& bobbinName) {
    auto requireSet = [&](bool isSet, const std::string& member) {
        if (!isSet) {
            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
                "The family '" + familyName + "' bobbin processor returned without setting " + member +
                " for bobbin '" + bobbinName + "'; a processed description is not assembled from an unset value.");
        }
    };
    requireSet(geometry.columnShape.has_value(), "columnShape");
    requireSet(geometry.columnThickness.has_value(), "columnThickness");
    requireSet(geometry.wallThickness.has_value(), "wallThickness");
    requireSet(geometry.columnWidth.has_value(), "columnWidth");
    requireSet(geometry.columnDepth.has_value(), "columnDepth");
    requireSet(geometry.windingWindow.has_value(), "the winding window");

    auto requirePositive = [&](double value, const std::string& member) {
        if (!(value > 0)) {
            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
                "Bobbin '" + bobbinName + "' (family '" + familyName + "') processed to a " + member + " of " +
                std::to_string(value) + " m; a central column has a positive extent.");
        }
    };
    requirePositive(geometry.columnWidth.value(), "columnWidth");
    requirePositive(geometry.columnDepth.value(), "columnDepth");

    CoreBobbinProcessedDescription processedDescription;
    processedDescription.set_column_shape(geometry.columnShape.value());
    processedDescription.set_column_thickness(geometry.columnThickness.value());
    processedDescription.set_wall_thickness(geometry.wallThickness.value());
    processedDescription.set_column_width(geometry.columnWidth.value());
    processedDescription.set_column_depth(geometry.columnDepth.value());
    processedDescription.get_mutable_winding_windows().push_back(geometry.windingWindow.value());
    processedDescription.set_coordinates(std::vector<double>({0, 0, 0}));
    return processedDescription;
}

}  // namespace OpenMagnetics
