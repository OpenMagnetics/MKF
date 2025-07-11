#pragma once

#include "json.hpp"

#include <MAS.hpp>
#include "Defaults.h"
#include "CorePiece.h"

#include <vector>

using namespace MAS;

using json = nlohmann::json;

namespace OpenMagnetics {

class Core : public MAS::MagneticCore {
  private:
    bool _includeMaterialData = false;

  public:
    Core(const json& j, bool includeMaterialData = false, bool includeProcessedDescription = true, bool includeGeometricalDescription = true);
    Core(const MagneticCore core);
    Core(const CoreShape shape, std::optional<CoreMaterial> material = std::nullopt);
    Core() = default;
    virtual ~Core() = default;

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
    void set_gap_length(double gapLength);
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
    
    double get_mass();

    double get_effective_length();
    double get_effective_area();
    double get_minimum_area();
    double get_effective_volume();

    std::vector<CoreGap> get_gapping() { return get_mutable_functional_description().get_gapping(); }
    double get_initial_permeability(double temperature = Defaults().ambientTemperature);
    static double get_initial_permeability(CoreMaterial coreMaterial, double temperature = Defaults().ambientTemperature);
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
    std::vector<VolumetricCoreLossesMethodType> get_available_core_losses_methods();
    bool can_be_used_for_filtering();
    CoreType get_type();
    bool fits(MaximumDimensions maximumDimensions, bool allowRotation=false);
    std::vector<double> get_maximum_dimensions();
    void set_material_initial_permeability(double value);
    void set_ground_gap(double gapLength);
    void set_distributed_gap(double gapLength, size_t numberGaps);
    void set_spacer_gap(double gapLength);
    void set_residual_gap();
};


void from_json(const json& j, std::vector<BhCycleElement>& v);
void to_json(json& j, const std::vector<BhCycleElement>& v);
void from_json(const json& j, std::vector<VolumetricLossesPoint>& v);
void to_json(json& j, const std::vector<VolumetricLossesPoint>& v);
void from_json(const json& j, std::vector<Permeability>& v);
void to_json(json& j, const std::vector<Permeability>& v);
void from_json(const json& j, std::vector<SteinmetzCoreLossesMethodRangeDatum>& v);
void to_json(json& j, const std::vector<SteinmetzCoreLossesMethodRangeDatum>& v);

inline void from_json(const json& j, std::vector<BhCycleElement>& v) {
    for (auto e : j) {
        BhCycleElement x(e);
        v.push_back(x);
    }
}

inline void to_json(json& j, const std::vector<BhCycleElement>& v) {
    j = json::array();
    for (auto x : v) {
        json e;
        to_json(e, x);
        j.push_back(e);
    }
}

inline void from_json(const json& j, std::vector<VolumetricLossesPoint>& v) {
    for (auto e : j) {
        VolumetricLossesPoint x(e);
        v.push_back(x);
    }
}

inline void to_json(json& j, const std::vector<VolumetricLossesPoint>& v) {
    j = json::array();
    for (auto x : v) {
        json e;
        to_json(e, x);
        j.push_back(e);
    }
}

inline void from_json(const json& j, std::vector<SteinmetzCoreLossesMethodRangeDatum>& v) {
    for (auto e : j) {
        SteinmetzCoreLossesMethodRangeDatum x(e);
        v.push_back(x);
    }
}

inline void to_json(json& j, const std::vector<SteinmetzCoreLossesMethodRangeDatum>& v) {
    j = json::array();
    for (auto x : v) {
        json e;
        to_json(e, x);
        j.push_back(e);
    }
}


} // namespace OpenMagnetics
