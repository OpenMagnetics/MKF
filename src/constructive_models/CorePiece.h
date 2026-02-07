#pragma once

#include "json.hpp"

#include <MAS.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <streambuf>
#include <vector>

using namespace MAS;

using json = nlohmann::json;

namespace OpenMagnetics {

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
};

void from_json(const json& j, OpenMagnetics::CorePiece& x);
void to_json(json& j, const OpenMagnetics::CorePiece& x);

} // namespace OpenMagnetics
