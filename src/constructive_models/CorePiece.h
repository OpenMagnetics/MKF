#pragma once

#include "json.hpp"

#include <MAS.hpp>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <optional>
#include <streambuf>
#include <vector>

using namespace MAS;

using json = nlohmann::json;

namespace OpenMagnetics {


/// Absolute volumes (m³) for each core part in the half-model.
struct CorePartVolumes {
    double centralColumn = 0.0;
    double lateralColumn = 0.0;   // one lateral column (symmetry)
    double topYoke        = 0.0;
    double bottomYoke     = 0.0;

    double total() const {
        return centralColumn + lateralColumn + topYoke + bottomYoke;
    }
};

/// Fractional loss distribution (sums to 1.0 for the modelled half-core).
struct CoreLossFractions {
    double centralColumn = 0.0;
    double lateralColumn = 0.0;
    double topYoke        = 0.0;
    double bottomYoke     = 0.0;
};

class CorePiece {
  private:
    std::vector<ColumnElement> columns;
    double depth;
    double height;
    double width;
    CoreShape shape;
    WindingWindowElement windingWindow;
    EffectiveParameters partialEffectiveParameters;

  public:
    virtual std::tuple<double, double, double> get_shape_constants() = 0;
    virtual std::tuple<double, double, double> get_shape_constants_iec63182();
    virtual void process_columns() = 0;
    virtual void process_winding_window() = 0;
    virtual void process_extra_data() = 0;

    // ABT #362: pieces whose magnetic circuit crosses TWO materials (e.g. a ferrite drum
    // closed by a magnetic-epoxy shell) expose their IEC shape constants split per material:
    // {c1_core, c2_core, c1_shell, c2_shell}. c1 alone drives the per-section reluctance
    // (mu applied section by section); c1 and c2 together give each material's own effective
    // parameters — le = c1^2/c2, Ae = c1/c2, Ve = c1^3/c2^2 — which the core-loss split needs
    // to price each material over ITS volume at ITS flux density. Single-material pieces
    // return nullopt (default).
    virtual std::optional<std::array<double, 4>> get_mixed_material_constants() { return std::nullopt; }

    // ABT #1002: a compression-moulded body can be pressed from up to THREE powders -- the post
    // the coil sits on, the cover moulded over the coil, and the base plate under it (the WE lists
    // of parts name them COR / COV / SUB, or Inner / Outer when base and cover are one pressing).
    // Pieces built that way expose their IEC 60205 sections grouped per REGION, in the order
    // functionalDescription.material lists the grades, so the inductance model can apply each
    // region's own permeability and the loss model can price each region's own volume. c1 and c2
    // are the same sums get_shape_constants() makes: adding the regions reproduces the piece.
    // Single-region pieces return nullopt (default).
    struct RegionShapeConstants {
        std::string name;
        double c1;
        double c2;
        double minimumArea;
    };
    virtual std::optional<std::vector<RegionShapeConstants>> get_region_shape_constants() { return std::nullopt; }

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
    const WindingWindowElement& get_winding_window() const { return windingWindow; }
    WindingWindowElement& get_mutable_winding_window() { return windingWindow; }
    void set_winding_window(const WindingWindowElement& value) { this->windingWindow = value; }

    const CoreShape get_shape() const { return shape; }
    CoreShape get_mutable_shape() { return shape; }
    void set_shape(CoreShape value) { this->shape = value; }

    const EffectiveParameters& get_partial_effective_parameters() const { return partialEffectiveParameters; }
    EffectiveParameters& get_mutable_partial_effective_parameters() { return partialEffectiveParameters; }
    void set_partial_effective_parameters(const EffectiveParameters& value) {
        this->partialEffectiveParameters = value;
    }

    static std::shared_ptr<CorePiece> factory(CoreShape shape, bool process=true);

    /**
     * @brief Whether factory() can build a piece for this shape family.
     *
     * ABT #307: MAS ships shape records for families whose geometry class does not
     * exist yet (UI, PQI, ...). The catalog loader uses this to leave those shapes
     * out instead of letting one unbuildable record abort an entire adviser sweep.
     * Single source of truth for factory()'s dispatch — keep the two in step.
     */
    static bool is_family_supported(CoreShapeFamily family);

    void process();
    
    // Thermal surface area calculations
    /**
     * @brief Get the central column's right face area (facing winding window)
     * @param columnIndex Index of the column (0 for central column)
     * @return Surface area in m²
     */
    double get_column_right_face_area(size_t columnIndex = 0);
    
    /**
     * @brief Get the central column's top face area
     * @param columnIndex Index of the column
     * @return Surface area in m²
     */
    double get_column_top_face_area(size_t columnIndex = 0);
    
    /**
     * @brief Get the central column's bottom face area
     * @param columnIndex Index of the column
     * @return Surface area in m²
     */
    double get_column_bottom_face_area(size_t columnIndex = 0);
    
    /**
     * @brief Get the yoke's interior face area (facing winding window)
     * @param isTopYoke True for top yoke, false for bottom yoke
     * @return Surface area in m²
     */
    double get_yoke_interior_face_area(bool isTopYoke);
    
    /**
     * @brief Get the yoke's exterior face area (facing away from winding)
     * @param isTopYoke True for top yoke, false for bottom yoke
     * @return Surface area in m²
     */
    double get_yoke_exterior_face_area(bool isTopYoke);
    
    /**
     * @brief Get the yoke's right face area (vertical face facing winding window)
     * @param isTopYoke True for top yoke, false for bottom yoke
     * @return Surface area in m²
     */
    double get_yoke_right_face_area(bool isTopYoke);
    
    /**
     * @brief Get the winding window height
     * @return Height in m
     */
    double get_winding_window_height();
    
    /**
     * @brief Get the winding window width
     * @return Width in m
     */
    double get_winding_window_width();
    
    /**
     * @brief Get the column width by index
     * @param columnIndex Index of the column
     * @return Width in m
     */
    double get_column_width(size_t columnIndex = 0);
    
    /**
     * @brief Get the column depth by index
     * @param columnIndex Index of the column
     * @return Depth in m
     */
    double get_column_depth(size_t columnIndex = 0);
    
    /**
     * @brief Get the column shape by index
     * @param columnIndex Index of the column
     * @return ColumnShape enum
     */
    ColumnShape get_column_shape(size_t columnIndex = 0);

    // ========================================================================
    // Volume-proportional core loss distribution
    // ========================================================================

    /**
     * @brief Calculate the cross-sectional area of a column, respecting its shape.
     *
     * Prefers the pre-computed ColumnElement::get_area() when > 0.
     * Falls back to a shape-aware formula (ROUND -> pi/4*w*d, OBLONG -> stadium,
     * RECTANGULAR/IRREGULAR -> w*d).
     *
     * @param columnIndex Index of the column (0 = central)
     * @return Cross-sectional area in m²
     */
    double calculate_column_cross_section(size_t columnIndex = 0);

    /**
     * @brief Compute half-model volumes (m³) for each core part.
     *
     * Uses a symmetry plane at the core centre (z = 0), so every
     * part stores HALF the physical depth.
     *
     * Column volume  = cross-section area × height / 2.
     * Yoke volume    = yokeWidth × yokeThickness × (coreDepth / 2),
     * where yokeWidth and yokeThickness match the dimensions used by
     * the thermal node initialisation in Temperature.cpp.
     *
     * @return CorePartVolumes with all fields populated
     */
    CorePartVolumes calculate_core_part_volumes();

    /**
     * @brief Normalise core part volumes into loss fractions that sum to 1.0.
     *
     * Assumes uniform volumetric loss density (valid for ungapped cores at
     * moderate flux densities).  For gapped cores the central-column fraction
     * may still be slightly overestimated, but the error is far smaller than
     * the previous hardcoded 40/20/10/10 split.
     *
     * @return CoreLossFractions (all fields in [0, 1], sum ≈ 1.0)
     */
    CoreLossFractions calculate_core_loss_fractions();
};

/**
 * The shape families CorePiece::factory can actually construct.
 *
 * This is the ENGINE's capability, and it is deliberately not the same question as
 * get_core_shape_families() in Utils.h, which reports the families that happen to
 * appear in the loaded shape database. A family with no catalogue shape is still
 * fully buildable from a custom shape, so a UI that offers families must ask this
 * one; asking the database instead silently hides every family nobody has published
 * a part for yet.
 *
 * @return the supported families, in declaration order
 */
std::vector<CoreShapeFamily> get_supported_core_shape_families();

/**
 * The dimensions a family's geometry reads and REQUIRES, from the CorePiece subclass that
 * reads them — not from whichever shapes are published, which answers nothing for a family
 * with no catalogue record (ABT #1007). Optional, guarded dimensions are deliberately absent;
 * get_shape_family_dimensions adds whatever the catalogue additionally carries.
 *
 * @throws std::runtime_error if the family has no declaration.
 */
std::vector<std::string> get_core_shape_family_required_dimensions(CoreShapeFamily family);

void from_json(const json& j, OpenMagnetics::CorePiece& x);
void to_json(json& j, const OpenMagnetics::CorePiece& x);

} // namespace OpenMagnetics
