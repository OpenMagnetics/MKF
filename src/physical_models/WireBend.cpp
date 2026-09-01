#include "physical_models/WireBend.h"

#include "Definitions.h"
#include "support/Exceptions.h"

#include <array>
#include <limits>
#include <cmath>
#include <numbers>

namespace OpenMagnetics {

namespace {

// Relative tolerance for the comparisons that pick a regime or a verdict. The quantities are
// millimetre-scale lengths held in metres, so an absolute floor keeps the comparison honest at
// small sizes without ever making a real 1 um difference invisible.
constexpr double kLengthTolerance = 1e-12;

struct MandrelRow {
    double conductorDiameter;  // [m]
    double mandrelDiameter;    // [m]
};

// IEC 60317-0-1:2013, Table 7 -- Heat shock. The wire is wound on this mandrel and THEN
// thermally shocked; the coating must not crack. Transcribed verbatim; the note under the table
// sends conductor diameters up to and including 0,140 mm to Table 6 instead, and requires an
// intermediate diameter to take the mandrel of the next larger tabulated size.
constexpr std::array<MandrelRow, 21> kRoundHeatShockMandrels = {{
    {0.160e-3, 0.250e-3}, {0.180e-3, 0.280e-3}, {0.200e-3, 0.315e-3}, {0.224e-3, 0.355e-3},
    {0.250e-3, 0.400e-3}, {0.280e-3, 0.630e-3}, {0.315e-3, 0.710e-3}, {0.355e-3, 0.800e-3},
    {0.400e-3, 0.900e-3}, {0.450e-3, 1.000e-3}, {0.500e-3, 1.120e-3}, {0.560e-3, 1.250e-3},
    {0.630e-3, 1.400e-3}, {0.710e-3, 1.600e-3}, {0.800e-3, 1.800e-3}, {0.900e-3, 2.000e-3},
    {1.000e-3, 2.240e-3}, {1.120e-3, 3.550e-3}, {1.250e-3, 4.000e-3}, {1.400e-3, 4.500e-3},
    {1.600e-3, 5.000e-3},
}};

// IEC 60317-0-1:2013, Table 6 -- Mandrel winding. Fine wire is pre-stretched and wound on a
// fixed 0,150 mm mandrel; from 0,140 mm upwards there is no pre-stretch and the mandrel IS the
// conductor. Above 1,600 mm the standard drops the mandrel test for a stretching test (8.2),
// so no bend radius exists to return.
constexpr double kRoundFixedMandrel = 0.150e-3;
constexpr double kRoundFixedMandrelMaximumDiameter = 0.140e-3;
constexpr double kRoundMandrelMaximumDiameter = 1.600e-3;

// IEC 60317-0-2, Table 6 -- Mandrel winding, rectangular copper. The mandrel is a multiple of
// the dimension that lies IN the bend plane, and §9 gives the heat-shock bend flatwise only.
constexpr double kRectangularMandrelFactor = 4.0;
constexpr double kRectangularWideMandrelFactor = 5.0;       // widths over 10 mm
constexpr double kRectangularWideThreshold = 10.0e-3;
constexpr double kRectangularHeatShockFlatwiseFactor = 6.0;

struct BendPlaneDimensions {
    double conducting;  // [m] bare conductor dimension lying in the bend plane
    double outer;       // [m] the same dimension over the insulation
};

double resolve_required(const std::optional<MAS::DimensionWithTolerance>& dimension,
                        const std::string& what) {
    if (!dimension) {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA,
                                   "WireBend: wire has no " + what +
                                       ", which the bend radius is defined against");
    }
    return resolve_dimensional_values(dimension.value());
}

BendPlaneDimensions bend_plane_dimensions(const Wire& wire, BendAxis axis) {
    if (wire.get_type() == WireType::ROUND) {
        Wire copy = wire;
        return {resolve_required(wire.get_conducting_diameter(), "conducting diameter"),
                Wire::calculate_outer_diameter(copy)};
    }
    if (wire.get_type() == WireType::RECTANGULAR) {
        Wire copy = wire;
        const double conductingWidth = resolve_required(wire.get_conducting_width(), "conducting width");
        const double conductingHeight = resolve_required(wire.get_conducting_height(), "conducting height");
        const double outerWidth = Wire::calculate_outer_width(copy);
        const double outerHeight = Wire::calculate_outer_height(copy);
        // The standard names the LARGER dimension "width" and the smaller one "thickness",
        // regardless of which axis MAS happens to store them on.
        const bool widthIsTheLargerDimension = conductingWidth >= conductingHeight;
        const bool wantTheLargerDimension = (axis == BendAxis::EDGEWISE);
        if (widthIsTheLargerDimension == wantTheLargerDimension) {
            return {conductingWidth, outerWidth};
        }
        return {conductingHeight, outerHeight};
    }
    throw InvalidInputException(
        ErrorCode::INVALID_WIRE_DATA,
        "WireBend: no bend requirement is standardised for this wire type; IEC 60317-0-1 and "
        "-0-2 cover round and rectangular copper only. Litz and foil need supplier data");
}

// The one case where the standards specify NO mandrel rather than an error: IEC 60317-0-2 gives
// the heat-shock bend flatwise only, so an edgewise heat-shock radius does not exist. Returning
// an empty optional keeps that absence distinguishable from a missing input.
std::optional<double> mandrel_diameter_if_specified(const Wire& wire, BendCriterion criterion,
                                                    BendAxis axis) {
    if (wire.get_type() == WireType::ROUND) {
        const double conductingDiameter =
            resolve_required(wire.get_conducting_diameter(), "conducting diameter");
        if (conductingDiameter <= kRoundFixedMandrelMaximumDiameter + kLengthTolerance) {
            // Table 7's own note: below 0,140 mm the heat-shock bend is Table 6's.
            return kRoundFixedMandrel;
        }
        if (conductingDiameter > kRoundMandrelMaximumDiameter + kLengthTolerance) {
            throw InvalidInputException(
                ErrorCode::INVALID_WIRE_DATA,
                "WireBend: IEC 60317-0-1 tabulates no mandrel above a 1,600 mm conductor "
                "(clause 8.2 replaces the winding test with a stretching test), so no minimum "
                "bend radius is defined for a conductor of " +
                    std::to_string(conductingDiameter) + " m");
        }
        if (criterion == BendCriterion::FLEXIBILITY) {
            return conductingDiameter;  // Table 6: the mandrel is the conductor itself
        }
        for (const auto& row : kRoundHeatShockMandrels) {
            // "for intermediate values, the mandrel of the next larger nominal diameter"
            if (conductingDiameter <= row.conductorDiameter + kLengthTolerance) {
                return row.mandrelDiameter;
            }
        }
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA,
                                    "WireBend: no heat-shock mandrel tabulated for a conductor of " +
                                        std::to_string(conductingDiameter) + " m");
    }

    if (wire.get_type() == WireType::RECTANGULAR) {
        const auto dimensions = bend_plane_dimensions(wire, axis);
        if (criterion == BendCriterion::FLEXIBILITY) {
            if (axis == BendAxis::EDGEWISE) {
                const double factor = dimensions.conducting > kRectangularWideThreshold
                                          ? kRectangularWideMandrelFactor
                                          : kRectangularMandrelFactor;
                return factor * dimensions.conducting;
            }
            return kRectangularMandrelFactor * dimensions.conducting;
        }
        if (axis == BendAxis::EDGEWISE) {
            return std::nullopt;  // §9 specifies the heat-shock bend flatwise only
        }
        return kRectangularHeatShockFlatwiseFactor * dimensions.conducting;
    }

    // Delegates the "not standardised" message, so there is one wording for it.
    bend_plane_dimensions(wire, axis);
    return std::nullopt;
}

BendVerdict verdict_for(double bendRadius, double flexibilityRadius,
                        const std::optional<double>& heatShockRadius) {
    if (bendRadius < flexibilityRadius - kLengthTolerance) {
        return BendVerdict::BELOW_FLEXIBILITY;
    }
    if (heatShockRadius && bendRadius < heatShockRadius.value() - kLengthTolerance) {
        return BendVerdict::BELOW_HEAT_SHOCK;
    }
    return BendVerdict::OK;
}

} // namespace

double WireBend::get_mandrel_diameter(const Wire& wire, BendCriterion criterion, BendAxis axis) {
    const auto mandrel = mandrel_diameter_if_specified(wire, criterion, axis);
    if (!mandrel) {
        throw InvalidInputException(
            ErrorCode::INVALID_WIRE_DATA,
            "WireBend: IEC 60317-0-2 specifies the heat-shock bend for rectangular wire flatwise "
            "only; there is no edgewise heat-shock mandrel to return");
    }
    return mandrel.value();
}

double WireBend::get_minimum_bend_radius(const Wire& wire, BendCriterion criterion, BendAxis axis) {
    const double mandrelDiameter = get_mandrel_diameter(wire, criterion, axis);
    const auto dimensions = bend_plane_dimensions(wire, axis);
    // IEC 60851-3 §5.1.1: the wire is wound ON the mandrel, so its centreline sits one mandrel
    // radius plus one wire radius from the axis.
    return 0.5 * mandrelDiameter + 0.5 * dimensions.outer;
}

WoundCorner WireBend::solve(double formerCornerRadius, double cornerHalfAngle, const Wire& wire,
                            double standoff, BendAxis axis, BendPolicy policy) {
    if (!std::isfinite(formerCornerRadius) || formerCornerRadius < 0) {
        throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
                                    "WireBend: former corner radius must be finite and not "
                                    "negative, got " +
                                        std::to_string(formerCornerRadius) + " m");
    }
    if (!std::isfinite(standoff) || standoff <= 0) {
        throw InvalidInputException(ErrorCode::CALCULATION_INVALID_INPUT,
                                    "WireBend: standoff must be finite and positive, got " +
                                        std::to_string(standoff) + " m");
    }
    if (!std::isfinite(cornerHalfAngle) || cornerHalfAngle <= 0 ||
        cornerHalfAngle > 0.5 * std::numbers::pi + kLengthTolerance) {
        throw InvalidInputException(ErrorCode::CALCULATION_INVALID_INPUT,
                                    "WireBend: corner half-angle must lie in (0, pi/2], got " +
                                        std::to_string(cornerHalfAngle) + " rad");
    }

    WoundCorner result{};
    result.flexibilityRadius = get_minimum_bend_radius(wire, BendCriterion::FLEXIBILITY, axis);
    result.heatShockRadius =
        mandrel_diameter_if_specified(wire, BendCriterion::HEAT_SHOCK, axis)
            ? std::optional<double>(get_minimum_bend_radius(wire, BendCriterion::HEAT_SHOCK, axis))
            : std::nullopt;

    double targetRadius = result.flexibilityRadius;
    if (policy == BendPolicy::CONSERVATIVE) {
        if (!result.heatShockRadius) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA,
                                       "WireBend: a conservative (heat-shock) placement was asked "
                                       "for, but this wire and bend axis have no heat-shock "
                                       "mandrel in the standard");
        }
        targetRadius = result.heatShockRadius.value();
    }
    result.mandrelDiameterUsed =
        get_mandrel_diameter(wire, policy == BendPolicy::CONSERVATIVE ? BendCriterion::HEAT_SHOCK
                                                                     : BendCriterion::FLEXIBILITY,
                             axis);

    const double available = formerCornerRadius + standoff;
    if (available >= targetRadius - kLengthTolerance) {
        // The wire follows the former: its corner is simply the former's corner offset by the
        // standoff, and the straight runs stay exactly where they were.
        result.regime = BendRegime::CONFORMING;
        result.bendRadius = available;
        result.faceStandoff = standoff;
        result.cornerCentreInset = 0.0;
    }
    else {
        // The wire cannot make that corner. It bends at its own limit, rides on the former's
        // corners and lifts off the flat faces. Tension pulls it in as far as the clearance to
        // the corner allows, which fixes the arc centre and hence the new standoff.
        result.regime = BendRegime::LIFTED;
        result.bendRadius = targetRadius;
        result.cornerCentreInset = targetRadius - formerCornerRadius - standoff;
        result.faceStandoff =
            standoff + (1.0 - std::sin(cornerHalfAngle)) * result.cornerCentreInset;
    }

    const auto dimensions = bend_plane_dimensions(wire, axis);
    result.outerFibreStrain = 0.5 * dimensions.conducting / result.bendRadius;
    result.verdict = verdict_for(result.bendRadius, result.flexibilityRadius, result.heatShockRadius);
    return result;
}

WoundCorner WireBend::evaluate(double bendRadius, const Wire& wire, BendAxis axis) {
    if (!std::isfinite(bendRadius) || bendRadius <= 0) {
        throw InvalidInputException(ErrorCode::CALCULATION_INVALID_INPUT,
                                    "WireBend: bend radius to evaluate must be finite and "
                                    "positive, got " +
                                        std::to_string(bendRadius) + " m");
    }
    WoundCorner result{};
    result.flexibilityRadius = get_minimum_bend_radius(wire, BendCriterion::FLEXIBILITY, axis);
    result.heatShockRadius =
        mandrel_diameter_if_specified(wire, BendCriterion::HEAT_SHOCK, axis)
            ? std::optional<double>(get_minimum_bend_radius(wire, BendCriterion::HEAT_SHOCK, axis))
            : std::nullopt;
    result.mandrelDiameterUsed = get_mandrel_diameter(wire, BendCriterion::FLEXIBILITY, axis);
    result.bendRadius = bendRadius;
    result.regime = BendRegime::CONFORMING;
    // An as-built radius carries no former, so there is no standoff or arc inset to report.
    // NaN rather than 0 so that reading them by mistake cannot pass for a real measurement.
    result.faceStandoff = std::numeric_limits<double>::quiet_NaN();
    result.cornerCentreInset = std::numeric_limits<double>::quiet_NaN();
    const auto dimensions = bend_plane_dimensions(wire, axis);
    result.outerFibreStrain = 0.5 * dimensions.conducting / bendRadius;
    result.verdict = verdict_for(bendRadius, result.flexibilityRadius, result.heatShockRadius);
    return result;
}

BendAxis WireBend::axis_from_bend_plane_dimension(const Wire& wire, bool bendPlaneIsWidth) {
    if (wire.get_type() == WireType::ROUND) {
        return BendAxis::ROUND;
    }
    if (wire.get_type() != WireType::RECTANGULAR) {
        // Same wording as everywhere else for an unsupported type.
        bend_plane_dimensions(wire, BendAxis::FLATWISE);
    }
    const double conductingWidth = resolve_required(wire.get_conducting_width(), "conducting width");
    const double conductingHeight = resolve_required(wire.get_conducting_height(), "conducting height");
    const double inBendPlane = bendPlaneIsWidth ? conductingWidth : conductingHeight;
    const double acrossBendPlane = bendPlaneIsWidth ? conductingHeight : conductingWidth;
    return inBendPlane >= acrossBendPlane ? BendAxis::EDGEWISE : BendAxis::FLATWISE;
}

std::string WireBend::to_string(BendRegime regime) {
    switch (regime) {
        case BendRegime::CONFORMING: return "CONFORMING";
        case BendRegime::LIFTED: return "LIFTED";
    }
    throw InvalidInputException("WireBend: unknown bend regime");
}

std::string WireBend::to_string(BendVerdict verdict) {
    switch (verdict) {
        case BendVerdict::OK: return "OK";
        case BendVerdict::BELOW_HEAT_SHOCK: return "BELOW_HEAT_SHOCK";
        case BendVerdict::BELOW_FLEXIBILITY: return "BELOW_FLEXIBILITY";
    }
    throw InvalidInputException("WireBend: unknown bend verdict");
}

} // namespace OpenMagnetics
