#include "physical_models/AxisymmetricWindow.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <numbers>
#include <string>
#include <tuple>

#include "support/Exceptions.h"

namespace OpenMagnetics {

namespace {

// Complete elliptic integrals K(m), E(m) of parameter m = k^2, Abramowitz & Stegun 17.3.34 and 17.3.36
// (|error| <= 2e-8), with m1 = 1 - m.
std::pair<double, double> complete_elliptic_integrals(double m) {
    const double m1 = std::max(1.0 - m, 1e-300);
    const double logarithm = -std::log(m1);
    const double K = (1.38629436112 + m1 * (0.09666344259 + m1 * (0.03590092383 + m1 * (0.03742563713 + m1 * 0.01451196212)))) +
                     (0.5 + m1 * (0.12498593597 + m1 * (0.06880248576 + m1 * (0.03328355346 + m1 * 0.00441787012)))) * logarithm;
    const double E = (1.0 + m1 * (0.44325141463 + m1 * (0.06260601220 + m1 * (0.04757383546 + m1 * 0.01736506451)))) +
                     m1 * (0.24998368310 + m1 * (0.09200180037 + m1 * (0.04069697526 + m1 * 0.00526449639))) * logarithm;
    return {K, E};
}

// Exponentially scaled modified Bessel functions, Abramowitz & Stegun 9.8.1-9.8.8 (relative error ~2e-7):
// I_n(x) e^-x and K_n(x) e^x, n = 0, 1, x > 0.
double scaled_i0(double x) {
    if (x <= 3.75) {
        const double t = (x / 3.75) * (x / 3.75);
        return std::exp(-x) * (1.0 + t * (3.5156229 + t * (3.0899424 + t * (1.2067492 + t * (0.2659732 + t * (0.0360768 + t * 0.0045813))))));
    }
    const double t = 3.75 / x;
    return (0.39894228 + t * (0.01328592 + t * (0.00225319 + t * (-0.00157565 + t * (0.00916281 + t * (-0.02057706 +
            t * (0.02635537 + t * (-0.01647633 + t * 0.00392377)))))))) / std::sqrt(x);
}

double scaled_i1(double x) {
    if (x <= 3.75) {
        const double t = (x / 3.75) * (x / 3.75);
        return std::exp(-x) * x * (0.5 + t * (0.87890594 + t * (0.51498869 + t * (0.15084934 + t * (0.02658733 + t * (0.00301532 + t * 0.00032411))))));
    }
    const double t = 3.75 / x;
    return (0.39894228 + t * (-0.03988024 + t * (-0.00362018 + t * (0.00163801 + t * (-0.01031555 + t * (0.02282967 +
            t * (-0.02895312 + t * (0.01787654 - t * 0.00420059)))))))) / std::sqrt(x);
}

double scaled_k0(double x) {
    if (x <= 2.0) {
        const double u = (x / 2) * (x / 2);
        const double i0 = scaled_i0(x) * std::exp(x);
        return std::exp(x) * (-std::log(x / 2) * i0 + (-0.57721566 + u * (0.42278420 + u * (0.23069756 + u * (0.03488590 + u * (0.00262698 +
                              u * (0.00010750 + u * 0.0000074)))))));
    }
    const double v = 2.0 / x;
    return (1.25331414 + v * (-0.07832358 + v * (0.02189568 + v * (-0.01062446 + v * (0.00587872 + v * (-0.00251540 + v * 0.00053208)))))) / std::sqrt(x);
}

double scaled_k1(double x) {
    if (x <= 2.0) {
        const double u = (x / 2) * (x / 2);
        const double i1 = scaled_i1(x) * std::exp(x);
        return std::exp(x) / x * (x * std::log(x / 2) * i1 + (1.0 + u * (0.15443144 + u * (-0.67278579 + u * (-0.18156897 + u * (-0.01919402 +
                                  u * (-0.00110404 - u * 0.00004686)))))));
    }
    const double v = 2.0 / x;
    return (1.25331414 + v * (0.23498619 + v * (-0.03655620 + v * (0.01504268 + v * (-0.00780353 + v * (0.00325614 - v * 0.00068245)))))) / std::sqrt(x);
}

// The pairs at one argument, sharing the asymptotic branch's 1 / sqrt(x) (the window's modes mostly take it).
std::pair<double, double> scaled_i0_i1(double x) {
    if (x <= 3.75) {
        return {scaled_i0(x), scaled_i1(x)};
    }
    const double t = 3.75 / x;
    const double inverseRoot = 1.0 / std::sqrt(x);
    return {(0.39894228 + t * (0.01328592 + t * (0.00225319 + t * (-0.00157565 + t * (0.00916281 + t * (-0.02057706 +
             t * (0.02635537 + t * (-0.01647633 + t * 0.00392377)))))))) * inverseRoot,
            (0.39894228 + t * (-0.03988024 + t * (-0.00362018 + t * (0.00163801 + t * (-0.01031555 + t * (0.02282967 +
             t * (-0.02895312 + t * (0.01787654 - t * 0.00420059)))))))) * inverseRoot};
}

std::pair<double, double> scaled_k0_k1(double x) {
    if (x <= 2.0) {
        return {scaled_k0(x), scaled_k1(x)};
    }
    const double v = 2.0 / x;
    const double inverseRoot = 1.0 / std::sqrt(x);
    return {(1.25331414 + v * (-0.07832358 + v * (0.02189568 + v * (-0.01062446 + v * (0.00587872 + v * (-0.00251540 + v * 0.00053208)))))) * inverseRoot,
            (1.25331414 + v * (0.23498619 + v * (-0.03655620 + v * (0.01504268 + v * (-0.00780353 + v * (0.00325614 - v * 0.00068245)))))) * inverseRoot};
}

// Modes kept: the correction at distance d from a cylinder decays as exp(-k d); 600 modes resolve d down to
// 12 B / (600 pi), about a fiftieth of the window height, and every point takes only the modes it needs.
constexpr size_t maximumModes = 600;
constexpr double decayCut = 14.0;  // exp(-14) ~ 8e-7

}  // namespace

bool AxisymmetricWindow::applies(const Core& core) {
    if (core.get_shape_family() == CoreShapeFamily::T) {
        return false;
    }
    const auto& processed = core.get_processed_description();
    if (!processed) {
        throw CoreNotProcessedException();
    }
    if (processed->get_winding_windows().size() != 1) {
        return false;
    }
    const auto& columns = processed->get_columns();
    const size_t mainColumnIndex = core.get_main_column_index();
    if (mainColumnIndex >= columns.size() || columns[mainColumnIndex].get_shape() != ColumnShape::ROUND) {
        return false;
    }
    size_t lateralColumns = 0;
    for (size_t columnIndex = 0; columnIndex < columns.size(); ++columnIndex) {
        if (columnIndex != mainColumnIndex && columns[columnIndex].get_type() == ColumnType::LATERAL) {
            ++lateralColumns;
        }
    }
    bool functionalGap = false;
    for (const auto& gap : core.get_functional_description().get_gapping()) {
        if (gap.get_type() == GapType::SUBTRACTIVE || gap.get_type() == GapType::ADDITIVE) {
            functionalGap = true;
        }
    }
    return functionalGap || lateralColumns <= 1;
}

bool AxisymmetricWindow::applies(const Core& core, const Coil& coil) {
    if (!applies(core)) {
        return false;
    }
    const auto& turns = coil.get_turns_description();
    if (!turns) {
        throw CoilNotProcessedException("Axisymmetric window: the coil has no turns description");
    }
    const auto window = from_core(core);
    const double tolerance = 1e-9 * (window.outer_radius() + window.top() - window.bottom());
    for (const auto& turn : turns.value()) {
        const auto& coordinates = turn.get_coordinates();
        if (coordinates.size() < 2 || !turn.get_dimensions() || turn.get_dimensions()->size() < 2) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Axisymmetric window: turn " + turn.get_name() + " has no 2D coordinates or dimensions");
        }
        const double width = turn.get_dimensions().value()[0];
        const double height = turn.get_dimensions().value()[1];
        const double angle = turn.get_rotation() ? turn.get_rotation().value() / 180 * std::numbers::pi : 0.0;
        const double halfX = 0.5 * (std::abs(width * std::cos(angle)) + std::abs(height * std::sin(angle)));
        const double halfY = 0.5 * (std::abs(width * std::sin(angle)) + std::abs(height * std::cos(angle)));
        const double radius = std::abs(coordinates[0]);
        if (radius - halfX < window.post_radius() - tolerance || radius + halfX > window.outer_radius() + tolerance ||
            coordinates[1] - halfY < window.bottom() - tolerance || coordinates[1] + halfY > window.top() + tolerance) {
            return false;
        }
    }
    return true;
}

AxisymmetricWindow AxisymmetricWindow::from_core(const Core& core) {
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException();
    }
    const auto& windingWindows = core.get_processed_description()->get_winding_windows();
    if (windingWindows.size() != 1 || !windingWindows[0].get_width() || !windingWindows[0].get_height()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Axisymmetric window: the core needs one winding window with width and height");
    }
    // The frame of CoilMesherCenterModel and of the neutral window's wall sheets: from the main column's face,
    // vertically centred.
    const double postRadius = core.get_columns()[core.get_main_column_index()].get_width() / 2;
    const double width = windingWindows[0].get_width().value();
    const double height = windingWindows[0].get_height().value();
    return AxisymmetricWindow(postRadius, postRadius + width, -height / 2, height / 2);
}

AxisymmetricWindow::AxisymmetricWindow(double postRadius, double outerRadius, double zBottom, double zTop)
    : _a(postRadius), _b(outerRadius), _zBottom(zBottom), _zTop(zTop), _height(zTop - zBottom), _modes(maximumModes) {
    if (!(postRadius > 0) || !(outerRadius > postRadius) || !(zTop > zBottom)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Axisymmetric window: needs 0 < post radius < outer radius and a positive height, got " +
                                    std::to_string(postRadius) + ", " + std::to_string(outerRadius) + ", " + std::to_string(zTop - zBottom));
    }
    _k.resize(_modes + 1);
    _postRatio.resize(_modes + 1);
    _outerRatio.resize(_modes + 1);
    _inverseScaledI0Outer.resize(_modes + 1);
    _inverseScaledK0Post.resize(_modes + 1);
    _scaledI0Post.resize(_modes + 1);
    _scaledK0Outer.resize(_modes + 1);
    for (size_t n = 1; n <= _modes; ++n) {
        const double k = static_cast<double>(n) * std::numbers::pi / _height;
        _k[n] = k;
        _postRatio[n] = scaled_i0(k * _a) / scaled_i0(k * _b) * std::exp(k * (_a - _b));
        _outerRatio[n] = scaled_k0(k * _b) / scaled_k0(k * _a) * std::exp(-k * (_b - _a));
        _inverseScaledI0Outer[n] = 1.0 / scaled_i0(k * _b);
        _inverseScaledK0Post[n] = 1.0 / scaled_k0(k * _a);
        _scaledI0Post[n] = scaled_i0(k * _a);
        _scaledK0Outer[n] = scaled_k0(k * _b);
    }
}

bool AxisymmetricWindow::contains(double r, double z) const {
    return r >= _a && r <= _b && z >= _zBottom && z <= _zTop;
}

std::pair<double, double> AxisymmetricWindow::ring_field(double current, double rs, double zs, double r, double z) {
    if (!(rs > 0) || !(r > 0)) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Axisymmetric window: a ring or a field point at r <= 0");
    }
    const double dz = z - zs;
    const double sumSquared = (rs + r) * (rs + r) + dz * dz;
    const double differenceSquared = (rs - r) * (rs - r) + dz * dz;
    if (!(differenceSquared > 0)) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Axisymmetric window: field point on a ring");
    }
    const auto [K, E] = complete_elliptic_integrals(4 * rs * r / sumSquared);
    const double prefactor = current / (2 * std::numbers::pi * std::sqrt(sumSquared));
    const double fieldZ = prefactor * (K + (rs * rs - r * r - dz * dz) / differenceSquared * E);
    const double fieldR = prefactor * dz / r * (-K + (rs * rs + r * r + dz * dz) / differenceSquared * E);
    return {fieldR, fieldZ};
}

std::pair<double, double> AxisymmetricWindow::ring_field_with_plate_images(double current, double rs, double zs, double r, double z) const {
    double fieldR = 0;
    double fieldZ = 0;
    // Image n of a point at zs: n B + (-1)^n (zs - zc) + zc, zc the window's mid-height.
    const double zCentre = 0.5 * (_zBottom + _zTop);
    for (int n = -planeImages; n <= planeImages; ++n) {
        const double zImage = n * _height + ((n % 2 == 0) ? 1.0 : -1.0) * (zs - zCentre) + zCentre;
        auto [ringR, ringZ] = ring_field(current, rs, zImage, r, z);
        fieldR += ringR;
        fieldZ += ringZ;
    }
    return {fieldR, fieldZ};
}

std::pair<double, double> AxisymmetricWindow::ring_minus_line_field(double current, double rs, double zs, double r, double z, double minimumDistance) {
    double dr = r - rs;
    double dz = z - zs;
    const double distance = std::hypot(dr, dz);
    if (distance < minimumDistance) {
        // Read the (smooth) difference on the circle of that radius, in the same direction.
        const double scale = minimumDistance / std::max(distance, 1e-300);
        if (distance > 0) {
            dr *= scale;
            dz *= scale;
        }
        else {
            dr = minimumDistance;
            dz = 0;
        }
    }
    auto [ringR, ringZ] = ring_field(current, rs, zs, rs + dr, zs + dz);
    // MKF's planar line current: Hx = -I (ys - y) / (2 pi rho^2), Hy = I (xs - x) / (2 pi rho^2).
    const double rhoSquared = dr * dr + dz * dz;
    const double lineR = current * dz / (2 * std::numbers::pi * rhoSquared);
    const double lineZ = -current * dr / (2 * std::numbers::pi * rhoSquared);
    return {ringR - lineR, ringZ - lineZ};
}

AxisymmetricWindow::WallMismatch AxisymmetricWindow::empty_mismatch() const {
    WallMismatch mismatch;
    mismatch.post.assign(_modes + 1, 0.0);
    mismatch.outer.assign(_modes + 1, 0.0);
    return mismatch;
}

void AxisymmetricWindow::add_ring(WallMismatch& mismatch, double current, double rs, double zs) const {
    if (!(rs > _a) || !(rs < _b)) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Axisymmetric window: a source at r = " + std::to_string(rs) +
                                   " outside the window (" + std::to_string(_a) + ", " + std::to_string(_b) + ")");
    }
    // Every image of the plates included: a ring stacked periodically is a solenoid, uniform H_z = I / B inside it.
    mismatch.postMean -= current / _height;
    const double phase = std::numbers::pi / _height * (zs - _zBottom);
    const double cosineStep = std::cos(phase);
    const double sineStep = std::sin(phase);
    double cosine = 1.0;
    double sine = 0.0;
    const double weight = 2 * current * rs / _height;
    // exp(k (a - rs)) and exp(k (rs - b)) for k = n pi / B, as running products.
    const double postStep = std::exp(std::numbers::pi / _height * (_a - rs));
    const double outerStep = std::exp(std::numbers::pi / _height * (rs - _b));
    double postDecay = 1.0;
    double outerDecay = 1.0;
    double largest = 0;
    for (size_t n = 1; n <= _modes; ++n) {
        const double nextCosine = cosine * cosineStep - sine * sineStep;
        sine = sine * cosineStep + cosine * sineStep;
        cosine = nextCosine;
        postDecay *= postStep;
        outerDecay *= outerStep;
        const double k = _k[n];
        // I0(ka) K1(k rs) and K0(kb) I1(k rs), scaled so that nothing overflows.
        const auto [i1, k1] = std::pair{scaled_i1(k * rs), scaled_k1(k * rs)};
        const double postValue = _scaledI0Post[n] * k1 * postDecay;
        const double outerValue = _scaledK0Outer[n] * i1 * outerDecay;
        mismatch.post[n] -= weight * k * postValue * cosine;
        mismatch.outer[n] += weight * k * outerValue * cosine;
        // The modes decay as exp(-k d), d the ring's distance to the nearer cylinder: past 1e-17 of the largest, nothing.
        const double magnitude = k * std::max(postValue, outerValue);
        largest = std::max(largest, magnitude);
        if (magnitude < 1e-17 * largest) {
            break;
        }
    }
}

void AxisymmetricWindow::add_wall_field(WallMismatch& mismatch, bool postFace, double zStart, double zEnd, double tangentialField) const {
    const double from = std::max(_zBottom, std::min(zStart, zEnd));
    const double to = std::min(_zTop, std::max(zStart, zEnd));
    if (!(to > from)) {
        return;
    }
    auto& coefficients = postFace ? mismatch.post : mismatch.outer;
    double& mean = postFace ? mismatch.postMean : mismatch.outerMean;
    mean += tangentialField * (to - from) / _height;
    for (size_t n = 1; n <= _modes; ++n) {
        const double k = _k[n];
        coefficients[n] += 2 * tangentialField / (_height * k) * (std::sin(k * (to - _zBottom)) - std::sin(k * (from - _zBottom)));
    }
}

AxisymmetricWindow::Correction AxisymmetricWindow::solve(const WallMismatch& mismatch) const {
    Correction correction;
    correction.p.assign(_modes + 1, 0.0);
    correction.q.assign(_modes + 1, 0.0);
    // The uniform mode (H_z constant, H_r = 0) moves both cylinders alike; by Ampere the two means it must meet
    // agree once every source of the window is in, so each part takes their average and the parts add up exactly.
    correction.uniform = 0.5 * (mismatch.postMean + mismatch.outerMean);
    for (size_t n = 1; n <= _modes; ++n) {
        // p I0(kr)/I0(kb) + q K0(kr)/K0(ka) meets the post coefficient at r = a and the outer one at r = b.
        const double determinant = _postRatio[n] * _outerRatio[n] - 1.0;
        correction.p[n] = (mismatch.post[n] * _outerRatio[n] - mismatch.outer[n]) / determinant;
        correction.q[n] = (mismatch.outer[n] * _postRatio[n] - mismatch.post[n]) / determinant;
    }
    // The modes that carry anything: past the last one above 1e-16 of the largest, every sum may stop.
    double largest = 0;
    for (size_t n = 1; n <= _modes; ++n) {
        largest = std::max(largest, std::abs(correction.p[n]) + std::abs(correction.q[n]));
    }
    correction.modes = 0;
    for (size_t n = _modes; n >= 1; --n) {
        if (std::abs(correction.p[n]) + std::abs(correction.q[n]) > 1e-16 * largest) {
            correction.modes = n;
            break;
        }
    }
    return correction;
}

AxisymmetricWindow::PointBasis AxisymmetricWindow::basis_at(double r, double z) const {
    PointBasis basis;
    // Each series to its own wall's decay: the I terms (the outer wall's reply) fall as exp(-k (b - r)), the K terms
    // (the post's) as exp(-k (r - a)). A point next to one wall needs many modes of that wall's series only.
    const auto modes_for = [&](double distance) {
        return distance > 0 ? std::min(_modes, static_cast<size_t>(std::ceil(decayCut * _height / (std::numbers::pi * distance)))) : _modes;
    };
    const size_t modesP = modes_for(std::max(_b - r, 0.0));
    const size_t modesQ = modes_for(std::max(r - _a, 0.0));
    const size_t modes = std::max(modesP, modesQ);
    basis.modes = modes;
    basis.radialP.assign(modes + 1, 0.0);
    basis.radialQ.assign(modes + 1, 0.0);
    basis.axialP.assign(modes + 1, 0.0);
    basis.axialQ.assign(modes + 1, 0.0);
    const double phase = std::numbers::pi / _height * (z - _zBottom);
    const double cosineStep = std::cos(phase);
    const double sineStep = std::sin(phase);
    // exp(k (r - b)) and exp(-k (r - a)) for k = n pi / B, as running products.
    const double growthStep = std::exp(std::numbers::pi / _height * (r - _b));
    const double decayStep = std::exp(-std::numbers::pi / _height * (r - _a));
    double growth = 1.0;
    double decay = 1.0;
    double cosine = 1.0;
    double sine = 0.0;
    for (size_t n = 1; n <= modes; ++n) {
        const double nextCosine = cosine * cosineStep - sine * sineStep;
        sine = sine * cosineStep + cosine * sineStep;
        cosine = nextCosine;
        growth *= growthStep;
        decay *= decayStep;
        const double x = _k[n] * r;
        // I_n(kr)/I0(kb) and K_n(kr)/K0(ka), scaled.
        if (n <= modesP) {
            const auto [i0, i1] = scaled_i0_i1(x);
            const double growing = growth * _inverseScaledI0Outer[n];
            basis.axialP[n] = i0 * growing * cosine;
            basis.radialP[n] = i1 * growing * sine;
        }
        if (n <= modesQ) {
            const auto [k0, k1] = scaled_k0_k1(x);
            const double decaying = decay * _inverseScaledK0Post[n];
            basis.axialQ[n] = k0 * decaying * cosine;
            basis.radialQ[n] = -k1 * decaying * sine;
        }
    }
    return basis;
}

std::pair<double, double> AxisymmetricWindow::correction_field(const Correction& correction, const PointBasis& basis) {
    double fieldR = 0;
    double fieldZ = correction.uniform;
    const size_t modes = std::min(basis.modes, correction.modes);
    for (size_t n = 1; n <= modes; ++n) {
        fieldR += correction.p[n] * basis.radialP[n] + correction.q[n] * basis.radialQ[n];
        fieldZ += correction.p[n] * basis.axialP[n] + correction.q[n] * basis.axialQ[n];
    }
    return {fieldR, fieldZ};
}

AxisymmetricWindow::WallParts AxisymmetricWindow::correction_field_by_wall(const Correction& correction, const PointBasis& basis) {
    WallParts parts;
    parts.post = {0.0, 0.5 * correction.uniform};
    parts.outer = {0.0, 0.5 * correction.uniform};
    const size_t modes = std::min(basis.modes, correction.modes);
    for (size_t n = 1; n <= modes; ++n) {
        parts.outer.first += correction.p[n] * basis.radialP[n];
        parts.outer.second += correction.p[n] * basis.axialP[n];
        parts.post.first += correction.q[n] * basis.radialQ[n];
        parts.post.second += correction.q[n] * basis.axialQ[n];
    }
    return parts;
}

std::pair<double, double> AxisymmetricWindow::correction_field(const Correction& correction, double r, double z) const {
    return correction_field(correction, basis_at(r, z));
}

AxisymmetricWindow::WallSources AxisymmetricWindow::wall_sources(const std::vector<WindowWallSheet>& sheetSources) const {
    WallSources sources;
    auto mismatch = empty_mismatch();
    const double tolerance = 1e-9 * (_b + _height);
    static constexpr double nodes[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
                                        0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
    static constexpr double weights[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
                                          0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
    for (const auto& sheet : sheetSources) {
        if (std::abs(sheet.start[0] - sheet.end[0]) < tolerance) {
            const double length = std::abs(sheet.end[1] - sheet.start[1]);
            if (!(length > 0)) {
                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Axisymmetric window: zero-length wall sheet");
            }
            const bool postFace = std::abs(sheet.start[0] - _a) < tolerance;
            if (!postFace && std::abs(sheet.start[0] - _b) >= tolerance) {
                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Axisymmetric window: a wall sheet at x = " + std::to_string(sheet.start[0]) +
                                           " is on neither cylinder of the window");
            }
            add_wall_field(mismatch, postFace, sheet.start[1], sheet.end[1], (postFace ? -1.0 : 1.0) * sheet.current / length);
        }
        else if (std::abs(sheet.start[1] - sheet.end[1]) < tolerance) {
            const double from = std::min(sheet.start[0], sheet.end[0]);
            const double to = std::max(sheet.start[0], sheet.end[0]);
            for (size_t node = 0; node < 8; ++node) {
                const double radius = 0.5 * (from + to) + 0.5 * (to - from) * nodes[node];
                const double current = sheet.current * 0.5 * weights[node];
                sources.plateRings.push_back({radius, sheet.start[1], current});
                add_ring(mismatch, current, radius, sheet.start[1]);
            }
        }
        else {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Axisymmetric window: a slanted wall sheet");
        }
    }
    sources.correction = solve(mismatch);
    return sources;
}

std::pair<double, double> AxisymmetricWindow::wall_field(const WallSources& sources, double pointX, double pointY) const {
    const double radius = std::abs(pointX);
    double fieldR = 0;
    double fieldZ = 0;
    for (const auto& ring : sources.plateRings) {
        auto [ringR, ringZ] = radius > 0 ? ring_field_with_plate_images(ring.current, ring.radius, ring.height, radius, pointY) : std::pair<double, double>{0.0, 0.0};
        fieldR += ringR;
        fieldZ += ringZ;
    }
    if (contains(radius, pointY)) {
        auto [correctionR, correctionZ] = correction_field(sources.correction, radius, pointY);
        fieldR += correctionR;
        fieldZ += correctionZ;
    }
    return {pointX < 0 ? -fieldR : fieldR, fieldZ};
}

std::vector<AxisymmetricWindow::WallFieldPart> AxisymmetricWindow::wall_field_parts(const WallSources& sources, double pointX, double pointY) const {
    const double radius = std::abs(pointX);
    const double side = pointX < 0 ? -1.0 : 1.0;
    std::vector<WallFieldPart> parts;
    parts.reserve(sources.plateRings.size() + 2);
    if (radius > 0) {
        for (const auto& ring : sources.plateRings) {
            auto [ringR, ringZ] = ring_field_with_plate_images(ring.current, ring.radius, ring.height, radius, pointY);
            parts.push_back({side * ring.radius, ring.height, side * ringR, ringZ});
        }
    }
    if (contains(radius, pointY)) {
        auto wallParts = correction_field_by_wall(sources.correction, basis_at(radius, pointY));
        parts.push_back({side * _a, pointY, side * wallParts.post.first, wallParts.post.second});
        parts.push_back({side * _b, pointY, side * wallParts.outer.first, wallParts.outer.second});
    }
    return parts;
}

double AxisymmetricWindow::ring_flux(double rs, double zs, double r, double z) {
    if (!(rs > 0)) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Axisymmetric window: a ring at r <= 0");
    }
    if (!(r > 0)) {
        return 0;
    }
    const double dz = z - zs;
    const double m = 4 * rs * r / ((rs + r) * (rs + r) + dz * dz);
    if (!(m < 1)) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Axisymmetric window: flux of a ring through itself");
    }
    const auto [K, E] = complete_elliptic_integrals(m);
    const double k = std::sqrt(m);
    constexpr double mu0 = 4e-7 * std::numbers::pi;
    return 2 * mu0 / k * std::sqrt(rs * r) * ((1 - m / 2) * K - E);
}

double AxisymmetricWindow::ring_flux_with_plate_images(double rs, double zs, double r, double z) const {
    double flux = 0;
    const double zCentre = 0.5 * (_zBottom + _zTop);
    for (int n = -planeImages; n <= planeImages; ++n) {
        const double zImage = n * _height + ((n % 2 == 0) ? 1.0 : -1.0) * (zs - zCentre) + zCentre;
        flux += ring_flux(rs, zImage, r, z);
    }
    return flux;
}

double AxisymmetricWindow::correction_flux(const Correction& correction, double r, double z) const {
    // psi with H_z = (1/r) d psi / dr and H_r = -(1/r) d psi / dz: c0 r^2 / 2 for the uniform mode and
    // (r / k) [p I1(kr) / I0(kb) - q K1(kr) / K0(ka)] cos(k (z - zBottom)) per mode.
    constexpr double mu0 = 4e-7 * std::numbers::pi;
    double psi = correction.uniform * r * r / 2;
    const double growthStep = std::exp(std::numbers::pi / _height * (r - _b));
    const double decayStep = std::exp(-std::numbers::pi / _height * (r - _a));
    double growth = 1.0;
    double decay = 1.0;
    const double phase = std::numbers::pi / _height * (z - _zBottom);
    const double cosineStep = std::cos(phase);
    const double sineStep = std::sin(phase);
    double cosine = 1.0;
    double sine = 0.0;
    for (size_t n = 1; n <= correction.modes; ++n) {
        const double nextCosine = cosine * cosineStep - sine * sineStep;
        sine = sine * cosineStep + cosine * sineStep;
        cosine = nextCosine;
        const double k = _k[n];
        growth *= growthStep;
        decay *= decayStep;
        const double growing = scaled_i1(k * r) * growth * _inverseScaledI0Outer[n];
        const double decaying = scaled_k1(k * r) * decay * _inverseScaledK0Post[n];
        const double term = r / k * (correction.p[n] * growing - correction.q[n] * decaying) * cosine;
        psi += term;
        if (std::abs(term) < 1e-16 * std::abs(psi) && n > 8) {
            break;
        }
    }
    return 2 * std::numbers::pi * mu0 * psi;
}

AxisymmetricWindow::FluxBasis AxisymmetricWindow::flux_basis_between(double r, double z, double referenceR, double referenceZ) const {
    constexpr double mu0 = 4e-7 * std::numbers::pi;
    const double scale = 2 * std::numbers::pi * mu0;
    FluxBasis basis;
    basis.uniform = scale * (r * r - referenceR * referenceR) / 2;
    basis.p.assign(_modes + 1, 0.0);
    basis.q.assign(_modes + 1, 0.0);
    // The same terms as correction_flux: (r / k) [p I1(kr) / I0(kb) - q K1(kr) / K0(ka)] cos(k (z - zBottom)).
    for (const auto& [radius, height, sign] : {std::tuple{r, z, 1.0}, std::tuple{referenceR, referenceZ, -1.0}}) {
        const double growthStep = std::exp(std::numbers::pi / _height * (radius - _b));
        const double decayStep = std::exp(-std::numbers::pi / _height * (radius - _a));
        const double phase = std::numbers::pi / _height * (height - _zBottom);
        const double cosineStep = std::cos(phase);
        const double sineStep = std::sin(phase);
        double growth = 1.0;
        double decay = 1.0;
        double cosine = 1.0;
        double sine = 0.0;
        for (size_t n = 1; n <= _modes; ++n) {
            const double nextCosine = cosine * cosineStep - sine * sineStep;
            sine = sine * cosineStep + cosine * sineStep;
            cosine = nextCosine;
            growth *= growthStep;
            decay *= decayStep;
            const double k = _k[n];
            const double factor = sign * scale * radius / k * cosine;
            basis.p[n] += factor * scaled_i1(k * radius) * growth * _inverseScaledI0Outer[n];
            basis.q[n] -= factor * scaled_k1(k * radius) * decay * _inverseScaledK0Post[n];
        }
    }
    return basis;
}

double AxisymmetricWindow::correction_flux(const Correction& correction, const FluxBasis& basis) {
    double flux = correction.uniform * basis.uniform;
    const size_t modes = std::min(basis.p.size() - 1, correction.modes);
    for (size_t n = 1; n <= modes; ++n) {
        flux += correction.p[n] * basis.p[n] + correction.q[n] * basis.q[n];
    }
    return flux;
}

double AxisymmetricWindow::wall_flux(const WallSources& sources, double r, double z) const {
    double flux = correction_flux(sources.correction, r, z);
    for (const auto& ring : sources.plateRings) {
        flux += ring.current * ring_flux_with_plate_images(ring.radius, ring.height, r, z);
    }
    return flux;
}

}  // namespace OpenMagnetics
