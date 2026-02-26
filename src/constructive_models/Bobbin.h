#pragma once

#include "json.hpp"

#include "constructive_models/Core.h"

#include <MAS.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <streambuf>
#include <vector>
#include "spline.h"
#include "support/Exceptions.h"

using json = nlohmann::json;

namespace OpenMagnetics {

inline tk::spline bobbinFillingFactorInterpWidth;
inline tk::spline bobbinFillingFactorInterpHeight;
inline tk::spline bobbinWindingWindowProportionInterpWidth;
inline tk::spline bobbinWindingWindowProportionInterpHeight;

inline double minBobbinWidth;
inline double maxBobbinWidth;
inline double minBobbinHeight;
inline double maxBobbinHeight;
inline double minWindingWindowWidth;
inline double maxWindingWindowWidth;
inline double minWindingWindowHeight;
inline double maxWindingWindowHeight;

class Bobbin : public MAS::Bobbin {
  private:
  public:
    Bobbin(const json& j, bool includeProcessedDescription = true) {
        from_json(j, *this);
        if (get_functional_description()) {
            if (includeProcessedDescription) {
                process_data();
            }
        }
        // else {
        //     if (!get_processed_description()) {
        //         throw std::runtime_error("Missing data in bobbin");
        //     }
        // }
    }

    Bobbin(const MAS::Bobbin bobbin) {
        set_functional_description(bobbin.get_functional_description());

        if (bobbin.get_processed_description()) {
            set_processed_description(bobbin.get_processed_description());
        }
        if (bobbin.get_distributors_info()) {
            set_distributors_info(bobbin.get_distributors_info());
        }
        if (bobbin.get_manufacturer_info()) {
            set_manufacturer_info(bobbin.get_manufacturer_info());
        }
        if (bobbin.get_name()) {
            set_name(bobbin.get_name());
        }
    }


    Bobbin() = default;
    virtual ~Bobbin() = default;
    void process_data();
    
    static double get_filling_factor(double windingWindowWidth, double windingWindowHeight);
    static std::vector<double> get_winding_window_dimensions(double coreWindingWindowWidth, double coreWindingWindowHeight);
    std::vector<double> get_winding_window_dimensions(size_t windingWindowIndex = 0);
    double get_winding_window_area(size_t windingWindowIndex = 0);
    WindingWindowShape get_winding_window_shape(size_t windingWindowIndex = 0);
    std::vector<double> get_winding_window_coordinates(size_t windingWindowIndex = 0);
    std::pair<double, double> get_column_and_wall_thickness(size_t windingWindowIndex = 0);
    WindingOrientation get_winding_window_sections_orientation(size_t windingWindowIndex = 0);
    CoilAlignment get_winding_window_sections_alignment(size_t windingWindowIndex = 0);
    static Bobbin create_quick_bobbin(double windingWindowHeight, double windingWindowWidth, ColumnShape shape = ColumnShape::ROUND);
    static Bobbin create_quick_bobbin(Core core, bool nullDimensions = false);
    static Bobbin create_quick_bobbin(Core core, double thickness);
    static Bobbin create_quick_bobbin(Core core, double wallThickness, double columnThickness);
    bool check_if_fits(double dimension, bool isHorizontalOrRadial = true, size_t windingWindowIndex = 0);
    void set_winding_orientation(WindingOrientation windingOrientation, size_t windingWindowIndex = 0);
    std::optional<WindingOrientation> get_winding_orientation(size_t windingWindowIndex = 0);
    std::vector<double> get_maximum_dimensions();
    
    // Thermal surface area calculations
    /**
     * @brief Get the column's right face area (facing winding window)
     * @param coreDepth The core depth (depth of the magnetic circuit)
     * @param windingWindowIndex Index of the winding window
     * @return Surface area in m²
     */
    double get_column_right_face_area(double coreDepth, size_t windingWindowIndex = 0);
    
    /**
     * @brief Get the column's top face area
     * @param coreDepth The core depth
     * @param windingWindowIndex Index of the winding window
     * @return Surface area in m²
     */
    double get_column_top_face_area(double coreDepth, size_t windingWindowIndex = 0);
    
    /**
     * @brief Get the column's bottom face area
     * @param coreDepth The core depth
     * @param windingWindowIndex Index of the winding window
     * @return Surface area in m²
     */
    double get_column_bottom_face_area(double coreDepth, size_t windingWindowIndex = 0);
    
    /**
     * @brief Get the yoke's interior face area (facing winding window, bottom of top yoke / top of bottom yoke)
     * @param coreDepth The core depth
     * @param isTopYoke True for top yoke, false for bottom yoke
     * @param windingWindowIndex Index of the winding window
     * @return Surface area in m²
     */
    double get_yoke_interior_face_area(double coreDepth, bool isTopYoke, size_t windingWindowIndex = 0);
    
    /**
     * @brief Get the yoke's exterior face area (facing away from winding, top of top yoke / bottom of bottom yoke)
     * @param coreDepth The core depth
     * @param isTopYoke True for top yoke, false for bottom yoke
     * @param windingWindowIndex Index of the winding window
     * @return Surface area in m²
     */
    double get_yoke_exterior_face_area(double coreDepth, bool isTopYoke, size_t windingWindowIndex = 0);
    
    /**
     * @brief Get the yoke's right face area (vertical face facing winding window)
     * @param wallThickness The bobbin wall thickness
     * @param coreDepth The core depth
     * @param windingWindowIndex Index of the winding window
     * @return Surface area in m²
     */
    double get_yoke_right_face_area(double wallThickness, double coreDepth, size_t windingWindowIndex = 0);
    
    /**
     * @brief Get the winding window height
     * @param windingWindowIndex Index of the winding window
     * @return Height in m
     */
    double get_winding_window_height(size_t windingWindowIndex = 0);
    
    /**
     * @brief Get the winding window width
     * @param windingWindowIndex Index of the winding window
     * @return Width in m
     */
    double get_winding_window_width(size_t windingWindowIndex = 0);
    
    /**
     * @brief Get the column width (full width including bobbin thickness)
     * @return Column width in m
     */
    double get_column_width();
    
    /**
     * @brief Get the column depth (full depth including bobbin thickness)
     * @return Column depth in m
     */
    double get_column_depth();
};



class BobbinDataProcessor{
    public:
        BobbinDataProcessor() = default;
        virtual ~BobbinDataProcessor() = default;
        virtual CoreBobbinProcessedDescription process_data(OpenMagnetics::Bobbin bobbin) = 0;

        static std::shared_ptr<BobbinDataProcessor> factory(OpenMagnetics::Bobbin bobbin);
};

} // namespace OpenMagnetics
