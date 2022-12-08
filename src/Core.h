#pragma once

#include <fstream>
#include <numbers>
#include <iostream>
#include <cmath>
#include <vector>
#include <filesystem>
#include <streambuf>
#include <nlohmann/json-schema.hpp>
#include <CoreTemplate.hpp>
#include <magic_enum.hpp>
#include "json.hpp"
using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;


namespace OpenMagnetics {
    using nlohmann::json;
    using ColumnShape = ShapeEnum;

    enum class DimensionalValues : int { MAXIMUM, NOMINAL, MINIMUM };

    template <OpenMagnetics::DimensionalValues preferredValue> 
    double resolve_dimensional_values(OpenMagnetics::Dimension dimensionValue);

    class CorePiece {
        private:
            std::vector<ColumnElement> columns;
            double depth;
            double height;
            double width;
            CoreShape shape;
            WindingWindowElement winding_window;
            EffectiveParameters partial_effective_parameters;
        protected:
            using dimensions_map = std::shared_ptr<std::map<std::string, Dimension>>;
        public:
            virtual std::tuple<double, double, double> get_shape_constants() = 0;
            virtual void process_columns() = 0;
            virtual void process_winding_window() = 0;
            virtual void process_extra_data() = 0;

            void flatten_dimensions()
            {
                for ( auto &dimension : *get_mutable_shape().get_dimensions() ) {
                    dimension.second = resolve_dimensional_values<OpenMagnetics::DimensionalValues::NOMINAL>(dimension.second);
                }
            }

            /**
             * List of columns in the piece
             */
            const std::vector<ColumnElement> & get_columns() const { return columns; }
            std::vector<ColumnElement> & get_mutable_columns() { return columns; }
            void set_columns(const std::vector<ColumnElement> & value) { this->columns = value; }

            /**
             * Total depth of the piece
             */
            const double & get_depth() const { return depth; }
            void set_depth(const double & value) { this->depth = value; }

            /**
             * Total height of the piece
             */
            const double & get_height() const { return height; }
            void set_height(const double & value) { this->height = value; }

            /**
             * Total width of the piece
             */
            const double & get_width() const { return width; }
            void set_width(const double & value) { this->width = value; }

            /**
             * List of winding windows, all elements in the list must be of the same type
             */
            const WindingWindowElement & get_winding_window() const { return winding_window; }
            WindingWindowElement & get_mutable_winding_window() { return winding_window; }
            void set_winding_window(const WindingWindowElement & value) { this->winding_window = value; }

            const CoreShape get_shape() const { return shape; }
            CoreShape get_mutable_shape() { return shape; }
            void set_shape(CoreShape value) { this->shape = value; }

            const EffectiveParameters & get_partial_effective_parameters() const { return partial_effective_parameters; }
            EffectiveParameters & get_mutable_partial_effective_parameters() { return partial_effective_parameters; }
            void set_partial_effective_parameters(const EffectiveParameters & value) { this->partial_effective_parameters = value; }


            static std::shared_ptr<CorePiece> factory(CoreShape shape);

            void process()
            {
                flatten_dimensions();
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

    void from_json(const json & j, OpenMagnetics::CorePiece & x);
    void to_json(json & j, const OpenMagnetics::CorePiece & x);


    class Core:public CoreTemplate {
        public:
        Core(const json & j) {
            from_json(j, *this);
            process_data();
            process_gap();
        }
        virtual ~Core() = default;

        std::shared_ptr<std::vector<GeometricalDescription>> create_geometrical_description();
        std::vector<ColumnElement> find_columns_by_type(ColumnType columnType);
        ColumnElement find_closest_column_by_coordinates(std::vector<double> coordinates);
        std::vector<CoreGap> find_gaps_by_type(GappingType gappingType);
        void scale_to_stacks(int64_t numberStacks);
        void process_gap();
        void process_data();
    };
}
