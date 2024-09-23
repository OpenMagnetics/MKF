#pragma once

#include "json.hpp"

#include <MAS.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>
using json = nlohmann::json;

namespace OpenMagnetics {

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

    virtual ~CorePiece() = default;

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

    static std::shared_ptr<CorePiece> factory(CoreShape shape, bool process=true);

    void process();
};

void from_json(const json& j, OpenMagnetics::CorePiece& x);
void to_json(json& j, const OpenMagnetics::CorePiece& x);

class CoreWrapper : public MagneticCore {
  private:
    bool _includeMaterialData = false;

  public:
    CoreWrapper(const json& j, bool includeMaterialData = false, bool includeProcessedDescription = true, bool includeGeometricalDescription = true);
    CoreWrapper(const MagneticCore core);
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
    bool is_gap_processed();
    bool process_gap();
    bool distribute_and_process_gap();
    void process_data();
    CoreMaterial resolve_material();
    static CoreMaterial resolve_material(CoreMaterialDataOrNameUnion coreMaterial);
    CoreShape resolve_shape();
    static CoreShape resolve_shape(CoreShapeDataOrNameUnion coreShape);
    std::vector<ColumnElement> get_columns();
    std::vector<WindingWindowElement> get_winding_windows();
    double get_depth();
    double get_height();
    double get_width();

    std::vector<CoreGap> get_gapping() { return get_mutable_functional_description().get_gapping(); }
    double get_initial_permeability(double temperature);
    static double get_initial_permeability(CoreMaterial coreMaterial, double temperature);
    double get_effective_permeability(double temperature);
    double get_reluctance(double temperature);
    double get_density();
    static double get_density(CoreMaterial coreMaterial);
    double get_resistivity(double temperature);
    static double get_resistivity(CoreMaterial coreMaterial, double temperature);
    double get_curie_temperature();
    static double get_curie_temperature(CoreMaterial coreMaterial);
    double get_remanence(double temperature);
    static double get_remanence(CoreMaterial coreMaterial, double temperature, bool returnZeroIfMissing = false);
    double get_coercive_force(double temperature);
    static double get_coercive_force(CoreMaterial coreMaterial, double temperature, bool returnZeroIfMissing = false);
    double get_magnetic_flux_density_saturation(bool proportion = true);
    static double get_magnetic_flux_density_saturation(CoreMaterial coreMaterial, bool proportion = true);
    double get_magnetic_flux_density_saturation(double temperature, bool proportion = true);
    static double get_magnetic_flux_density_saturation(CoreMaterial coreMaterial, double temperature, bool proportion = true);
    double get_magnetic_field_strength_saturation(double temperature);
    static double get_magnetic_field_strength_saturation(CoreMaterial coreMaterial, double temperature);
    CoreShapeFamily get_shape_family();
    std::string get_shape_name();
    std::string get_material_name();
    int64_t get_number_stacks();
    std::vector<CoreLossesMethodType> get_available_core_losses_methods();
    OpenMagnetics::CoreType get_type();
    bool fits(MaximumDimensions maximumDimensions, bool allowRotation=false);
};
} // namespace OpenMagnetics
