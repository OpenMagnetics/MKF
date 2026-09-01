#pragma once

#include "constructive_models/Wire.h"

#include <optional>
#include <string>

namespace OpenMagnetics {

// MINIMUM BEND RADIUS OF A WOUND CONDUCTOR, AND WHAT HAPPENS WHEN THE FORMER IS TIGHTER.
//
// A wire wound onto a former bends around its corners exactly the way it bends around the
// mandrel of the standard flexibility test, so the standards specify the former's corner
// radius directly. IEC 60851-3 §5.1.1 winds the wire on a polished mandrel of the diameter
// given by the requirement standard, which makes the wire's CENTRELINE radius
//
//      R = mandrelDiameter / 2 + outerDiameter / 2
//
// and, since a turn laid on a former stands off by outerDiameter / 2 from its surface, the
// outerDiameter / 2 cancels: the former's corner radius must be at least the mandrel's radius.
//
// Requirement tables, all transcribed from the standards themselves:
//   * round copper      IEC 60317-0-1:2013 Table 6 (mandrel winding) and Table 7 (heat shock)
//   * rectangular copper IEC 60317-0-2 Table 6 (mandrel winding) and §9 (heat shock)
//
// Two criteria, because the standards give two and they differ by more than 2x:
//   FLEXIBILITY is the tightest bend the insulation is required to survive AT ALL (for round
//               copper from 0,140 mm to 1,600 mm the mandrel is the conductor diameter itself);
//   HEAT_SHOCK  is the bend it must survive and THEN be thermally shocked -- the criterion that
//               matches a real coil, which gets soldered, varnish-baked and cycled after winding.
//
// When the former's corner is tighter than the wire can bend, the wire does not follow it. It
// rides on the corners and lifts off the flat faces, which is the second regime below.
enum class BendCriterion : int {
    FLEXIBILITY,  // IEC 60317-0-1 Table 6 / IEC 60317-0-2 Table 6
    HEAT_SHOCK    // IEC 60317-0-1 Table 7 / IEC 60317-0-2 §9
};

// Which of the conductor's dimensions lies IN the bend plane. A round wire is rotationally
// symmetric so it has only one case; a rectangular one bends far more easily on its thickness
// (flatwise) than on its width (edgewise), and the two mandrels differ by the aspect ratio.
enum class BendAxis : int {
    ROUND,
    FLATWISE,  // bent on the thickness: the SMALL dimension lies in the bend plane
    EDGEWISE   // bent on the width:     the LARGE dimension lies in the bend plane
};

// Which criterion the GEOMETRY is resolved against. A winding machine pulls the wire against
// the former, so what physically limits the bend is the wire's absolute limit, not the
// comfortable one -- STANDARD_LIMIT is therefore the default. CONSERVATIVE places the winding
// as a manufacturer targeting heat-shock survival would. The verdict is reported against BOTH
// criteria either way, so a design that only works at the limit says so.
enum class BendPolicy : int {
    STANDARD_LIMIT,
    CONSERVATIVE
};

enum class BendRegime : int {
    CONFORMING,  // the wire follows the former's corner; the turn sits where it always did
    LIFTED       // the wire bridges: it rides the corners and stands off the flat faces
};

enum class BendVerdict : int {
    OK,                 // at or above the heat-shock bend radius
    BELOW_HEAT_SHOCK,   // legal to wind, but the insulation is not qualified to be heated after
    BELOW_FLEXIBILITY   // tighter than the wire is qualified to bend at all -- not windable
};

struct WoundCorner {
    BendRegime regime;
    BendVerdict verdict;
    double bendRadius;          // centreline radius actually achieved [m]
    double faceStandoff;        // offset of the straight runs from the former's face [m]
    double cornerCentreInset;   // arc centre inset from the former's corner centre [m]
    double outerFibreStrain;    // conductor half-dimension / bendRadius, dimensionless
    double flexibilityRadius;   // minimum centreline radius, FLEXIBILITY criterion [m]
    std::optional<double> heatShockRadius;  // ditto HEAT_SHOCK; absent when the standard
                                            // specifies no heat-shock mandrel for this case
    double mandrelDiameterUsed; // the mandrel the geometry was resolved against [m]
};

class WireBend {
    public:
        // Mandrel diameter [m] from the requirement standard, applying its own rule that an
        // intermediate conductor size takes the mandrel of the NEXT LARGER tabulated size.
        // Throws for wire types the standards do not cover (litz, foil, planar) and for
        // conductor sizes outside the tables -- never substitutes a value.
        static double get_mandrel_diameter(const Wire& wire, BendCriterion criterion, BendAxis axis);

        // Minimum centreline bend radius [m] = mandrel / 2 + outer dimension in the bend plane / 2.
        static double get_minimum_bend_radius(const Wire& wire, BendCriterion criterion, BendAxis axis);

        // THE SOLVER. Given the former's corner and the wire, return the geometry that actually
        // results. cornerHalfAngle is half the interior angle between the two straight runs the
        // corner joins: pi/4 at the corner of a rectangular column, and the lift-off term
        // vanishes as it approaches pi/2 (no corner at all).
        //
        //   available = formerCornerRadius + standoff
        //   available >= Rmin -> CONFORMING: bend = available, standoff unchanged
        //   available <  Rmin -> LIFTED:     bend = Rmin,
        //                                    standoff += (1 - sin(halfAngle)) * (Rmin - Rc - s)
        //
        // The lifted solution is the one tension produces: the wire is pulled in as far as it
        // can go, which leaves it touching the former's corners and clear of its flat faces.
        static WoundCorner solve(double formerCornerRadius,
                                 double cornerHalfAngle,
                                 const Wire& wire,
                                 double standoff,
                                 BendAxis axis,
                                 BendPolicy policy = BendPolicy::STANDARD_LIMIT);

        // Judge a bend radius that already exists (an as-built geometry) against the standards,
        // without proposing a different one. This is what audits current geometry.
        static WoundCorner evaluate(double bendRadius, const Wire& wire, BendAxis axis);

        // The bend axis implied by which of a rectangular wire's dimensions lies in the bend
        // plane. Round wires always answer ROUND.
        static BendAxis axis_from_bend_plane_dimension(const Wire& wire, bool bendPlaneIsWidth);

        static std::string to_string(BendRegime regime);
        static std::string to_string(BendVerdict verdict);
};

} // namespace OpenMagnetics
