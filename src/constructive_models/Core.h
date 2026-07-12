#pragma once

#include "json.hpp"
#include "constructive_models/MasMigration.h"

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

    // Internal cache used by the const overloads of resolve_material() /
    // resolve_shape(). The non-const overloads continue to memoise into
    // get_mutable_functional_description() (preserving the legacy
    // side-effect that many callers rely on, e.g.
    // `std::get<CoreShape>(get_functional_description().get_shape())`).
    // The const overloads populate this cache only — they do not mutate
    // MAS state — which is what lets read-only adviser code accept a
    // `const Core&` and still resolve material/shape names cheaply.
    // Mutators that change the underlying material/shape must invalidate.
    mutable std::optional<CoreMaterial> _cachedResolvedMaterial;
    mutable std::optional<CoreShape>    _cachedResolvedShape;

  public:
    Core(json j, bool includeMaterialData = false, bool includeProcessedDescription = true, bool includeGeometricalDescription = true);
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
    bool is_gapping_misaligned();
    bool is_gap_processed();
    bool process_gap();
    bool distribute_and_process_gap();
    void set_gap_length(double gapLength);
    void process_data();
    CoreMaterial resolve_material();
    CoreMaterial resolve_material() const;
    static CoreMaterial resolve_material(CoreMaterialDataOrNameUnion coreMaterial);
    CoreShape resolve_shape();
    CoreShape resolve_shape() const;
    static CoreShape resolve_shape(CoreShapeDataOrNameUnion coreShape);
    std::vector<ColumnElement> get_columns() const;
    std::vector<WindingWindowElement> get_winding_windows() const;
    WindingWindowElement get_winding_window(size_t windingWindowIndex = 0) const;
    // Index of the main column: the first CENTRAL column, or the first column if there
    // is no central one (e.g. U/UT cores). This is the column the historical single
    // winding window wraps.
    size_t get_main_column_index() const;
    // Index of the column the turns placed in the given winding window are wound
    // around: the window's `column` edge when present, else the main column (the
    // schema-documented default). Throws if the window index or a stamped column
    // index is out of range.
    size_t get_winding_window_column_index(size_t windingWindowIndex = 0) const;
    double get_depth() const;
    double get_height() const;
    double get_width() const;
    
    double get_mass();

    double get_effective_length() const;
    double get_effective_area() const;
    double get_minimum_area() const;
    double get_effective_volume() const;
    std::string get_reference() const;

    const std::vector<CoreGap>& get_gapping() const { return get_functional_description().get_gapping(); }
    // Coating thickness per surface, in m. 0 when the core has no coating. An explicit
    // thickness (CoreCoating object) is used as-is; a name-only coating (legacy string
    // form) resolves to the datasheet default thickness for that coating type.
    double get_coating_thickness() const;
    // Relative permittivity of the coating dielectric, from its insulation material.
    // Explicit CoreCoating.material wins; otherwise the datasheet default material for
    // the coating type. Throws if the core has no coating or the type has no material.
    double get_coating_relative_permittivity() const;
    // For an uncoated toroid, whether the default coating type is parylene (small cores)
    // rather than epoxy (larger cores), selected by the outer diameter.
    bool get_default_toroid_coating_is_parylene() const;
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
    CoreShapeFamily get_shape_family() const;
    std::string get_material_family() const;
    std::string get_shape_name() const;
    std::string get_material_name() const;
    MAS::MagneticApplication resolve_material_application();
    static MAS::MagneticApplication resolve_material_application(CoreMaterial& coreMaterial);
    MAS::MagneticApplication guess_material_application();
    static MAS::MagneticApplication guess_material_application(CoreMaterial coreMaterial);
    static MAS::MagneticApplication guess_material_application(std::string coreMaterialName);
    bool check_material_application(MAS::MagneticApplication application);
    static bool check_material_application(CoreMaterial coreMaterial, MAS::MagneticApplication application);
    int64_t get_number_stacks() const;
    std::vector<VolumetricCoreLossesMethodType> get_available_core_losses_methods() const;
    static std::vector<VolumetricCoreLossesMethodType> get_available_core_losses_methods(CoreMaterial coreMaterial);
    CoreType get_type() const;
    bool fits(MaximumDimensions maximumDimensions, bool allowRotation=false);
    std::vector<double> get_maximum_dimensions();
    void set_type(CoreType coreType);
    void set_number_stacks(int64_t numberStacks);
    void set_material(CoreMaterial coreMaterial);
    void set_shape(CoreShape coreShape);
    void set_gapping(std::vector<CoreGap> coreGapping);
    void set_material_initial_permeability(double value);

    void set_ground_gapping(double gapLength);
    void set_distributed_gapping(double gapLength, size_t numberGaps);
    void set_spacer_gapping(double gapLength);
    void set_residual_gapping();

    std::vector<CoreGap> create_ground_gapping(double gapLength);
    static std::vector<CoreGap> create_ground_gapping(double gapLength, size_t numberColumns);
    std::vector<CoreGap> create_distributed_gapping(double gapLength, size_t numberGaps);
    static std::vector<CoreGap> create_distributed_gapping(double gapLength, size_t numberGaps, size_t numberColumns);
    std::vector<CoreGap> create_spacer_gapping(double gapLength);
    static std::vector<CoreGap> create_spacer_gapping(double gapLength, size_t numberColumns);
    std::vector<CoreGap> create_residual_gapping();
    static std::vector<CoreGap> create_residual_gapping(size_t numberColumns);

    static Core create_quick_core(std::string coreShapeName, std::string coreMaterialName, std::vector<CoreGap> gapping = {}, int64_t numberStacks = 1);

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
