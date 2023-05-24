#pragma once

#include "Utils.h"
#include "json.hpp"

#include <MAS.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <nlohmann/json-schema.hpp>
#include <numbers>
#include <streambuf>
#include <vector>
using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;

namespace OpenMagnetics {
using nlohmann::json;

class CorePiece {
  private:
    std::vector<ColumnElement> columns;
    double depth;
    double height;
    double width;
    CoreShape shape;
    WindingWindowElement winding_window;
    EffectiveParameters partial_effective_parameters;

  public:
    virtual std::tuple<double, double, double> get_shape_constants() = 0;
    virtual void process_columns() = 0;
    virtual void process_winding_window() = 0;
    virtual void process_extra_data() = 0;

    /**
     * List of columns in the piece
     */
    const std::vector<ColumnElement>& get_columns() const { return columns; }
    std::vector<ColumnElement>& get_mutable_columns() { return columns; }
    void set_columns(const std::vector<ColumnElement>& value) { this->columns = value; }

    /**
     * Total depth of the piece
     */
    const double& get_depth() const { return depth; }
    void set_depth(const double& value) { this->depth = value; }

    /**
     * Total height of the piece
     */
    const double& get_height() const { return height; }
    void set_height(const double& value) { this->height = value; }

    /**
     * Total width of the piece
     */
    const double& get_width() const { return width; }
    void set_width(const double& value) { this->width = value; }

    /**
     * List of winding windows, all elements in the list must be of the same type
     */
    const WindingWindowElement& get_winding_window() const { return winding_window; }
    WindingWindowElement& get_mutable_winding_window() { return winding_window; }
    void set_winding_window(const WindingWindowElement& value) { this->winding_window = value; }

    const CoreShape get_shape() const { return shape; }
    CoreShape get_mutable_shape() { return shape; }
    void set_shape(CoreShape value) { this->shape = value; }

    const EffectiveParameters& get_partial_effective_parameters() const { return partial_effective_parameters; }
    EffectiveParameters& get_mutable_partial_effective_parameters() { return partial_effective_parameters; }
    void set_partial_effective_parameters(const EffectiveParameters& value) {
        this->partial_effective_parameters = value;
    }

    static std::shared_ptr<CorePiece> factory(CoreShape shape);

    void process() {
        set_shape(flatten_dimensions(get_mutable_shape()));
        process_winding_window();
        process_columns();
        process_extra_data();

        auto [c1, c2, minimumArea] = get_shape_constants();
        json pieceEffectiveParameters;
        pieceEffectiveParameters["effectiveLength"] = pow(c1, 2) / c2;
        pieceEffectiveParameters["effectiveArea"] = c1 / c2;
        pieceEffectiveParameters["effectiveVolume"] = pow(c1, 3) / pow(c2, 2);
        pieceEffectiveParameters["minimumArea"] = minimumArea;
        set_partial_effective_parameters(pieceEffectiveParameters);
    }
};

void from_json(const json& j, OpenMagnetics::CorePiece& x);
void to_json(json& j, const OpenMagnetics::CorePiece& x);

class CoreWrapper : public MagneticCore {
  private:
    bool _includeMaterialData = false;

  public:
    CoreWrapper(const json& j, bool includeMaterialData = false, bool includeProcessedDescription = true, bool includeGeometricalDescription = true) {
        _includeMaterialData = includeMaterialData;
        from_json(j, *this);
        
        if (includeProcessedDescription) {
            process_data();
            process_gap();
        }

        if (!get_geometrical_description() && includeGeometricalDescription) {
            auto geometricalDescription = create_geometrical_description();

            set_geometrical_description(geometricalDescription);
        }
    }
    CoreWrapper(const MagneticCore core) {
        set_functional_description(core.get_functional_description());

        if (core.get_geometrical_description()) {
            set_geometrical_description(core.get_geometrical_description());
        }
        if (core.get_processed_description()) {
            set_processed_description(core.get_processed_description());
        }
        if (core.get_distributors_info()) {
            set_distributors_info(core.get_distributors_info());
        }
        if (core.get_manufacturer_info()) {
            set_manufacturer_info(core.get_manufacturer_info());
        }
        if (core.get_name()) {
            set_name(core.get_name());
        }
    }
    CoreWrapper() = default;
    virtual ~CoreWrapper() = default;

    std::optional<std::vector<CoreGeometricalDescriptionElement>> create_geometrical_description();
    std::vector<ColumnElement> find_columns_by_type(ColumnType columnType);
    ColumnElement find_closest_column_by_coordinates(std::vector<double> coordinates);
    int find_closest_column_index_by_coordinates(std::vector<double> coordinates);
    int find_exact_column_index_by_coordinates(std::vector<double> coordinates);
    std::vector<CoreGap> find_gaps_by_type(GapType gappingType);
    std::vector<CoreGap> find_gaps_by_column(ColumnElement column);
    void scale_to_stacks(int64_t numberStacks);
    bool is_gapping_missaligned();
    void process_gap();
    void distribute_and_process_gap();
    void process_data();
    CoreMaterial get_material();
    std::vector<ColumnElement> get_columns();
    std::vector<WindingWindowElement> get_winding_windows();

    std::vector<CoreGap> get_gapping() { return get_mutable_functional_description().get_gapping(); }
    double get_magnetic_flux_density_saturation(double temperature, bool proportion = true);
    double get_magnetic_fielda_strength_saturation(double temperature);
    double get_remanence(double temperature);
    double get_coercive_force(double temperature);
};
} // namespace OpenMagnetics
