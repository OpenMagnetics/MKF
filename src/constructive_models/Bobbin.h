#pragma once

#include "json.hpp"
#include "constructive_models/MasMigration.h"

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

// ABT #113: per-thread memo caches, lazily rebuilt per thread from the
// (frozen) bobbin catalog — lock-free and semantically transparent. The
// min/max trackers below are built together with the splines and must stay
// in the same thread_local group.
inline thread_local tk::spline bobbinFillingFactorInterpWidth;
inline thread_local tk::spline bobbinFillingFactorInterpHeight;
inline thread_local tk::spline bobbinWindingWindowProportionInterpWidth;
inline thread_local tk::spline bobbinWindingWindowProportionInterpHeight;

inline thread_local double minBobbinWidth;
inline thread_local double maxBobbinWidth;
inline thread_local double minBobbinHeight;
inline thread_local double maxBobbinHeight;
inline thread_local double minWindingWindowWidth;
inline thread_local double maxWindingWindowWidth;
inline thread_local double minWindingWindowHeight;
inline thread_local double maxWindingWindowHeight;

// Smallest wall/column thickness observed in the bobbin database. Used as a hard
// lower bound when constructing quick bobbins for cores outside the interpolator's
// training range, where the proportion clamp would otherwise yield ~µm walls.
inline thread_local double minBobbinWallThickness;
inline thread_local double minBobbinColumnThickness;

class Bobbin : public MAS::Bobbin {
  private:
  public:
    Bobbin(json j, bool includeProcessedDescription = true) {
        OpenMagnetics::compat::migrate_pre_1_0(j);
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
    /**
     * @brief Expand a catalogue pinout into the individual pins it describes (WP2, ABT #1171).
     *
     * A `pinout` is a footprint recipe - "14 pins, 6 and 8 to a row, rows 10.16 mm apart,
     * pitch 3.81 and 2.54, 5.08 between the middle pair" - and nothing downstream can draw,
     * export or simulate a recipe. This turns it into `MAS::Pin`s carrying real coordinates,
     * referred to the centre of the main column, exactly as `bobbin.json` `$defs/pin` defines.
     *
     * FRAME (MVB++ concentric frame, proposal section 0.7): the column axis is Y, the winding
     * window is radial in X and deep in Z.
     *   VERTICAL  (column axis normal to the board): the pins leave the pin rail below the
     *             bottom flange along -Y.
     *             A row runs along X; the rows sit at z = -+ rowDistance/2.
     *   HORIZONTAL(column axis parallel to the board): the part lies on its side and the pins
     *             leave the rails of the two END flanges downwards, along -Z. A row runs along X; the rows
     *             sit at y = -+ rowDistance/2, i.e. one row per flange, which is why a
     *             horizontal bobbin's row distance scales with the core's window HEIGHT
     *             (ETD 44 -> 35.8 mm, ETD 49 -> 40.4 mm) and a vertical one's with its depth.
     * `rowDistance` is the FULL row-to-row distance, so a row sits at half of it - the reading
     * that matches every record in the database and Shulin's "Row Pitch" column alike.
     *
     * NUMBERING: "1".."N", counter-clockwise seen from the pin side, starting at the far end of
     * row 0 (the row at the negative side). Row 0 is numbered along +X, row 1 back along -X, so
     * the numbers walk the footprint as a ring - the de-facto convention on every DIP former.
     *
     * THROWS rather than inventing: a pinout that states a pin COUNT and no pitch or row
     * distance describes no geometry at all, and half a footprint is worse than none. Same for
     * a missing `pinDescription`: `pin.dimensions` is required by the schema, and a pin needs
     * its diameter and length to be a solid.
     *
     * WHERE THE PINS START (ABT #1207): at the outer face of the bobbin's PIN RAIL (the
     * standoff), never at a flange face. A flange face lies INSIDE the core window, so a pin
     * hung from it runs straight through the core's back plate; the rail is the part of the
     * former that reaches past the core, and only the record's own dimensions say where it is.
     * `get_pin_rail_distance` reads that from the record, and throws when it cannot.
     *
     * @param pinout           the catalogue pinout
     * @param orientation      mounting orientation; it alone decides which way the pins leave
     * @param pinRailDistance  distance from the centre of the main column to the pin rail's
     *                         outer face, measured along the direction the pins leave (-Y for
     *                         vertical, -Z for horizontal); the pin starts there
     * @return the pins, in name order "1".."N"
     *
     * The proposal's draft signature also took `columnWidth`; nothing in the placement reads it
     * (a row's extent comes from its own pitches), and an unused parameter is a lie about what
     * the function depends on, so it is not taken.
     */
    static std::vector<MAS::Pin> expand_pinout(const MAS::Pinout& pinout,
                                               MAS::OrientationEnum orientation,
                                               double pinRailDistance);

    /**
     * @brief Distance from the centre of the main column to the outer face of the pin rail
     *        (the standoff the pins leave from), along the pin direction (ABT #1207).
     *
     * Read ONLY from dimensions whose meaning has been checked against the vendor drawing:
     *   PQ, VERTICAL: `c - H1/2`. On every Miles-Platts PQ former drawing the catalogue's PQ
     *     records come from (PQ0010..PQ0080, e.g. PQ0040 = PQ 26/25: c = 1.033 in, H1 = 0.609 in)
     *     `c` is the overall height from the TOP flange's outer face to the pin standoff, and
     *     `H1` the flange-to-flange height; the bobbin is centred on the column, so the top
     *     flange face is at +H1/2 and the standoff at H1/2 - c.
     * Every other family/orientation has no dimension in MAS that locates the rail (the
     * letters `a`/`b`/`c` mean different things per family and vendor - `c` on an E bobbin is
     * not an overall height), so this throws naming what is missing rather than guessing.
     *
     * @throws InvalidInputException when the record has no orientation, lacks a label the rule
     *         needs, or belongs to a family/orientation with no verified rule.
     */
    static double get_pin_rail_distance(const MAS::BobbinFunctionalDescription& functionalDescription);

    /**
     * @brief The processed pin with this name, e.g. "7".
     * @throws when the bobbin has not been processed, carries no pins, or has no such pin.
     */
    MAS::Pin get_pin(const std::string& name);

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

    /**
     * @brief Radius of the corners of the central column, which is what a wire bends around
     *        while being wound onto it (see WireBend: the former's corner plays exactly the
     *        role of the mandrel in the IEC 60317 flexibility test).
     *
     * Derived where the shape already defines it -- a round column's corner radius is its own
     * radius, an oblong one's is half its depth -- and read from the MAS datum for rectangular
     * and irregular columns. A synthesised bobbin that carries no datum falls back to the
     * injection-moulding rule (inside radius >= 0.5 x wall thickness), which is the same kind of
     * modelled-but-sourced value create_quick_bobbin already uses for the wall thickness itself.
     * @return Column corner radius in m
     */
    double get_column_corner_radius();

    /**
     * @brief Half the interior angle between the two straight runs that a turn's corner joins:
     *        pi/4 on a rectangular column, pi/2 where there is no corner at all.
     * @return Corner half-angle in radians
     */
    double get_column_corner_half_angle();
};



class BobbinDataProcessor{
    public:
        BobbinDataProcessor() = default;
        virtual ~BobbinDataProcessor() = default;
        virtual CoreBobbinProcessedDescription process_data(OpenMagnetics::Bobbin bobbin) = 0;

        static std::shared_ptr<BobbinDataProcessor> factory(OpenMagnetics::Bobbin bobbin);
};

} // namespace OpenMagnetics
