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
};

void from_json(const json& j, OpenMagnetics::CorePiece& x);
void to_json(json& j, const OpenMagnetics::CorePiece& x);

} // namespace OpenMagnetics
