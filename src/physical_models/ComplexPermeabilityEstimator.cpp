#include "physical_models/ComplexPermeabilityEstimator.h"
#include "support/Logger.h"
#include "support/Utils.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <map>
#include <numbers>
#include <magic_enum.hpp>
#include <sstream>

// ABT #1861: THE APPROVED ESTIMATE of complex permeability outside a material's tabulated range
// (Alf, 2026-10-14). See ComplexPermeabilityEstimator.h for the method and its validation. This file
// is a port of the proposal's validated Python (.cperm-estimator/cperm.py, holdout.py: rank_donors,
// fit_donor, predict_per_point, extend_tails, tail_slope, hybrid); keep the two in step.

namespace OpenMagnetics {

namespace {

// Numerical settings of the method (none is a physical constant; every physical quantity is
// learned from the catalogue at call time):
// - the frequency shift tried per donor: +-1.5 decades in 0.025-decade steps;
constexpr double kShiftRangeDecades = 1.5;
constexpr double kShiftStepDecades = 0.025;
constexpr int kNumberShifts = 121;  // 2 * 1.5 / 0.025 + 1
// - the target's own data matched (Alf, 2026-10-14: each side uses the variant that did best on its
//   own hold-out): ABOVE the data, all its points on the known side of the gap edge, each weighted
//   exp(-d / 1 decade) with d its distance from the edge in decades (the neighbours give the
//   tendency, the material's own data drive the weight; hold-out: 0.5, 1 and 2 decades tried, 1
//   kept); BELOW the data, its points within 1 decade of the edge, uniformly. Points of the other part
//   inside the gap weigh 1;
constexpr double kMatchWeightDecades = 1.0;
constexpr double kMatchWindowDecades = 1.0;
// - mu' joins its own edge value through an additive offset that fades linearly in log-frequency to
//   0 over 1 decade from the edge (Alf, 2026-10-14);
constexpr double kEdgeOffsetFadeDecades = 1.0;
// - the join is made with the donors as they are just inside the gap, 1e-6 decades past the edge: a
//   donor whose shifted table ends exactly at the edge counts at the edge but not beyond it, and the
//   join at the edge point itself would leave a jump right after it (3C97's mu' at 30 MHz: 0.47);
constexpr double kEdgeJoinStepDecades = 1e-6;
// - mu' points inside the gap join the match only from 100 Hz up (the hold-out's sweep starts at
//   1 kHz; a decade below it is still matched);
constexpr double kLowestMatchedRealFrequency = 100;
// - a donor must reach at least 80 % (and at least 2) of a channel's points at a shift;
constexpr double kMinimumMatchedFraction = 0.8;
constexpr size_t kMinimumMatchedPoints = 2;
// - the estimate at a point is the median of the first 5 donors (in match order) that reach it,
//   among the 15 best, and needs at least 3 of them;
constexpr size_t kNeighbours = 5;
constexpr size_t kDonorsConsidered = 15;
constexpr size_t kMinimumNeighbours = 3;
// - the Snoek screen: a member whose mu_i * f(mu'' peak) lies more than 2 decades from its class
//   median is not a donor (the class needs 5 members with a peak to have a median);
constexpr double kSnoekScreenDecades = 2.0;
constexpr size_t kSnoekScreenMinimumMembers = 5;
// - a mu'' peak counts only strictly inside the data (2 points from each end, 20 % above both ends);
constexpr double kPeakProminence = 1.2;
// - near-duplicates of the target (the same curve under two names, e.g. ACME D37 / K151) are not
//   donors: initial permeability within 0.002 decades and mu'' within 0.005 decades (median over
//   more than 5 of the target's points);
constexpr double kDuplicateInitialPermeabilityDecades = 0.002;
constexpr double kDuplicateCurveDecades = 0.005;
constexpr size_t kDuplicateMinimumPoints = 5;
// - tail slopes: the donor's outermost decade, at least 3 points spanning more than a factor 3;
constexpr size_t kTailMinimumPoints = 3;
constexpr double kTailMinimumSpan = 3;
// - log-frequency tolerance at a table's ends.
constexpr double kLogTolerance = 1e-9;
// - anchored Kramers-Kronig mu' (ABT #1870; port of .cperm-estimator/physical/kk.py and knn_kk.py
//   im_curve): the own mu'' on a 200-per-decade log grid, padded 1 decade each side with the slope of the
//   power law fitted to the table's outermost decade (the high slope capped at -0.5: the tail must fall),
//   starting from the table's end value, grid points within 1e-12 f^2 of the singularity left out of the
//   quadrature;
constexpr double kKramersKronigGridPerDecade = 200;
constexpr double kKramersKronigPadDecades = 1;
constexpr double kKramersKronigEdgeFitDecades = 1;
constexpr double kKramersKronigMaximumHighSlope = -0.5;
constexpr double kKramersKronigSingularTolerance = 1e-12;
// - KK hands over to the neighbours over the last 0.25 decade of the own mu'' table, linearly in log f;
//   beyond the table the neighbours alone (Alf, 2026-10-14; options study: the 1-decade offset fade carried
//   KK's end level into the neighbours and flipped the sign of 5 records' mu');
constexpr double kKramersKronigBlendDecades = 0.25;
// - KK only for FERRITE NI_ZN (Alf, 2026-10-14). For MnZn the anchor offset C = mu'(f_a) - KK(f_a) is
//   mu'(infinity) - 1, and on TDK MnZn records whose mu' and mu'' curves were digitised a few hundredths of
//   a decade apart (N41, N45, N51, N92, N96, T36, T57, T66) C sets an implausible high-frequency mu' (N51:
//   +465 at 1 GHz) that no check on the overlap can see; MnZn keeps the neighbours (ABT #1861).
constexpr MaterialComposition kKramersKronigComposition = MaterialComposition::NI_ZN;

struct CatalogueMaterial {
    std::string name;
    std::string manufacturer;
    double initialPermeability;
    ComplexPermeabilityTable real;
    ComplexPermeabilityTable imaginary;
    std::optional<double> lossPeakFrequency;
};

double median(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    const size_t n = values.size();
    return n % 2 == 1 ? values[n / 2] : (values[n / 2 - 1] + values[n / 2]) / 2.0;
}

// Least-squares straight line y = slope * x + intercept (numpy polyfit degree 1).
std::pair<double, double> linear_fit(const std::vector<double>& xs, const std::vector<double>& ys) {
    double meanX = 0, meanY = 0;
    for (size_t index = 0; index < xs.size(); ++index) {
        meanX += xs[index];
        meanY += ys[index];
    }
    meanX /= static_cast<double>(xs.size());
    meanY /= static_cast<double>(xs.size());
    double sxy = 0, sxx = 0;
    for (size_t index = 0; index < xs.size(); ++index) {
        sxy += (xs[index] - meanX) * (ys[index] - meanY);
        sxx += (xs[index] - meanX) * (xs[index] - meanX);
    }
    double slope = sxy / sxx;
    return {slope, meanY - slope * meanX};
}

// Linear interpolation, clamped to the end values outside x (numpy interp). x strictly increasing.
double interpolate_clamped(const std::vector<double>& x, const std::vector<double>& y, double at) {
    if (at <= x.front()) {
        return y.front();
    }
    if (at >= x.back()) {
        return y.back();
    }
    size_t index = static_cast<size_t>(std::upper_bound(x.begin(), x.end(), at) - x.begin());
    double fraction = (at - x[index - 1]) / (x[index] - x[index - 1]);
    return y[index - 1] + fraction * (y[index] - y[index - 1]);
}

std::string format_frequency(double frequency) {
    std::ostringstream text;
    text << std::setprecision(3);
    if (frequency >= 1e9) {
        text << frequency / 1e9 << " GHz";
    }
    else if (frequency >= 1e6) {
        text << frequency / 1e6 << " MHz";
    }
    else if (frequency >= 1e3) {
        text << frequency / 1e3 << " kHz";
    }
    else {
        text << frequency << " Hz";
    }
    return text.str();
}

std::string part_symbol(ComplexPermeabilityPart part) {
    return part == ComplexPermeabilityPart::REAL ? "mu'" : "mu''";
}

std::vector<PermeabilityPoint> tabulated_points(const MAS::Permeability& part) {
    if (std::holds_alternative<std::vector<PermeabilityPoint>>(part)) {
        return std::get<std::vector<PermeabilityPoint>>(part);
    }
    return {std::get<PermeabilityPoint>(part)};
}

// Sorted by frequency, exact duplicates averaged; mu'' keeps only its positive points. A frequency
// that is not a finite positive number is malformed data and throws.
ComplexPermeabilityTable clean_table(const std::vector<PermeabilityPoint>& points, bool positiveOnly, const std::string& materialName) {
    std::map<double, std::pair<double, size_t>> accumulated;
    for (const auto& point : points) {
        if (!point.get_frequency()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Complex permeability point without frequency in " + materialName);
        }
        double frequency = point.get_frequency().value();
        if (!std::isfinite(frequency) || frequency <= 0 || !std::isfinite(point.get_value())) {
            throw MaterialException(ErrorCode::MATERIAL_INVALID_PROPERTY,
                "Core material '" + materialName + "' has a malformed complex permeability point in MAS (frequency " +
                std::to_string(frequency) + " Hz, value " + std::to_string(point.get_value()) + "). Fix the record in MAS.");
        }
        if (positiveOnly && point.get_value() <= 0) {
            continue;
        }
        auto& entry = accumulated[frequency];
        entry.first += point.get_value();
        entry.second += 1;
    }
    ComplexPermeabilityTable table;
    for (const auto& [frequency, sumAndCount] : accumulated) {
        table.logFrequencies.push_back(std::log10(frequency));
        table.values.push_back(sumAndCount.first / static_cast<double>(sumAndCount.second));
        if (positiveOnly) {
            table.logValues.push_back(std::log10(table.values.back()));
        }
    }
    return table;
}

// Interpolation linear in log f (mu': linear in value; mu'': linear in log value). Disengaged
// outside the table (with a 1e-9-decade tolerance at its ends).
std::optional<double> interpolate(const ComplexPermeabilityTable& table, double logFrequency, bool logValues) {
    const auto& x = table.logFrequencies;
    if (x.size() < 2 || logFrequency < x.front() - kLogTolerance || logFrequency > x.back() + kLogTolerance) {
        return std::nullopt;
    }
    const auto& y = logValues ? table.logValues : table.values;
    double clamped = std::clamp(logFrequency, x.front(), x.back());
    auto upper = std::upper_bound(x.begin(), x.end(), clamped);
    size_t index = upper == x.end() ? x.size() - 1 : static_cast<size_t>(upper - x.begin());
    if (index == 0) {
        index = 1;
    }
    double fraction = (clamped - x[index - 1]) / (x[index] - x[index - 1]);
    double value = y[index - 1] + fraction * (y[index] - y[index - 1]);
    return logValues ? std::pow(10.0, value) : value;
}

// The material's initial permeability as the catalogue states it at 25 C: its points without a
// frequency (else those at its lowest frequency), the one closest to 25 C (no temperature = 25 C).
double catalogue_initial_permeability(const CoreMaterial& material) {
    std::vector<PermeabilityPoint> points;
    const auto& initial = material.get_permeability().get_initial();
    if (std::holds_alternative<PermeabilityPoint>(initial)) {
        points.push_back(std::get<PermeabilityPoint>(initial));
    }
    else {
        points = std::get<std::vector<PermeabilityPoint>>(initial);
    }
    if (points.empty()) {
        throw MaterialDataMissingException(material.get_name(), "initial permeability");
    }
    std::vector<PermeabilityPoint> chosen;
    for (const auto& point : points) {
        if (!point.get_frequency() || point.get_frequency().value() == 0) {
            chosen.push_back(point);
        }
    }
    if (chosen.empty()) {
        double lowest = std::numeric_limits<double>::max();
        for (const auto& point : points) {
            lowest = std::min(lowest, point.get_frequency().value());
        }
        for (const auto& point : points) {
            if (point.get_frequency().value() == lowest) {
                chosen.push_back(point);
            }
        }
    }
    constexpr double referenceTemperature = 25;
    const PermeabilityPoint* best = nullptr;
    double bestDistance = std::numeric_limits<double>::max();
    for (const auto& point : chosen) {
        double distance = std::fabs(point.get_temperature().value_or(referenceTemperature) - referenceTemperature);
        if (distance < bestDistance) {
            bestDistance = distance;
            best = &point;
        }
    }
    return best->get_value();
}

// Frequency of the mu'' maximum, only when it lies strictly inside the data.
std::optional<double> loss_peak_frequency(const ComplexPermeabilityTable& imaginary) {
    const auto& values = imaginary.values;
    if (values.size() < 5) {
        return std::nullopt;
    }
    size_t peak = static_cast<size_t>(std::max_element(values.begin(), values.end()) - values.begin());
    if (peak < 2 || peak > values.size() - 3) {
        return std::nullopt;
    }
    if (values[peak] < kPeakProminence * std::max(values.front(), values.back())) {
        return std::nullopt;
    }
    return std::pow(10.0, imaginary.logFrequencies[peak]);
}

std::optional<CatalogueMaterial> catalogue_material(const CoreMaterial& material) {
    const auto& complex = material.get_permeability().get_complex();
    if (!complex) {
        return std::nullopt;
    }
    CatalogueMaterial entry;
    entry.name = material.get_name();
    entry.manufacturer = material.get_manufacturer_info().get_name();
    entry.real = clean_table(tabulated_points(complex->get_real()), false, entry.name);
    entry.imaginary = clean_table(tabulated_points(complex->get_imaginary()), true, entry.name);
    if (entry.real.values.size() < 2 || entry.imaginary.values.size() < 2) {
        return std::nullopt;
    }
    entry.initialPermeability = catalogue_initial_permeability(material);
    entry.lossPeakFrequency = loss_peak_frequency(entry.imaginary);
    return entry;
}

bool is_near_duplicate(const CatalogueMaterial& target, const CatalogueMaterial& donor) {
    if (!(std::fabs(std::log10(donor.initialPermeability / target.initialPermeability)) < kDuplicateInitialPermeabilityDecades)) {
        return false;
    }
    std::vector<double> differences;
    for (size_t index = 0; index < target.imaginary.values.size(); ++index) {
        auto donorValue = interpolate(donor.imaginary, target.imaginary.logFrequencies[index], true);
        if (donorValue) {
            differences.push_back(std::fabs(std::log10(donorValue.value() / target.imaginary.values[index])));
        }
    }
    return differences.size() > kDuplicateMinimumPoints && median(differences) < kDuplicateCurveDecades;
}

// A matching channel: the target's own points (log10 f, log10 value). "magnitude" is mu' (its
// positive points); "imaginary" is mu''.
struct Channel {
    std::vector<double> logFrequencies;
    std::vector<double> logValues;
    std::vector<double> weights;
};

struct Channels {
    Channel magnitude;
    Channel imaginary;
};

bool beyond_edge(double logFrequency, double logEdge, ComplexPermeabilityGapSide side) {
    return side == ComplexPermeabilityGapSide::BELOW ? logFrequency < logEdge : logFrequency > logEdge;
}

// The weight of one of the target's own points on the known side of the edge (0: not matched).
// Gap above the data: exp(-d / kMatchWeightDecades). Gap below it: 1 within kMatchWindowDecades.
double known_side_weight(double logFrequency, double logEdge, ComplexPermeabilityGapSide side) {
    double distance = side == ComplexPermeabilityGapSide::BELOW ? logFrequency - logEdge : logEdge - logFrequency;
    if (distance < -kLogTolerance) {
        return 0;
    }
    if (side == ComplexPermeabilityGapSide::BELOW) {
        return distance <= kMatchWindowDecades + kLogTolerance ? 1 : 0;
    }
    return std::exp(-std::max(distance, 0.0) / kMatchWeightDecades);
}

// The target's data used to fit each donor: its points on the known side of the edge, weighted
// as known_side_weight says, plus the
// other part's points inside the gap when that part is tabulated there (realMissing / imaginaryMissing
// say which parts the target lacks beyond the edge).
Channels match_channels(const CatalogueMaterial& target, double logEdge, ComplexPermeabilityGapSide side,
                        bool realMissing, bool imaginaryMissing) {
    Channels channels;
    const double logLowestMatchedReal = std::log10(kLowestMatchedRealFrequency);
    for (size_t index = 0; index < target.real.values.size(); ++index) {
        double logFrequency = target.real.logFrequencies[index];
        double value = target.real.values[index];
        double weight = known_side_weight(logFrequency, logEdge, side);
        if (!realMissing && beyond_edge(logFrequency, logEdge, side) && logFrequency >= logLowestMatchedReal) {
            weight = 1;
        }
        if (weight > 0 && value > 0) {
            channels.magnitude.logFrequencies.push_back(logFrequency);
            channels.magnitude.logValues.push_back(std::log10(value));
            channels.magnitude.weights.push_back(weight);
        }
    }
    for (size_t index = 0; index < target.imaginary.values.size(); ++index) {
        double logFrequency = target.imaginary.logFrequencies[index];
        double weight = known_side_weight(logFrequency, logEdge, side);
        if (!imaginaryMissing && beyond_edge(logFrequency, logEdge, side)) {
            weight = 1;
        }
        if (weight > 0) {
            channels.imaginary.logFrequencies.push_back(logFrequency);
            channels.imaginary.logValues.push_back(target.imaginary.logValues[index]);
            channels.imaginary.weights.push_back(weight);
        }
    }
    return channels;
}

struct DonorMatch {
    double residual;
    double shiftDecades;
    std::optional<double> magnitudeOffset;  // log10 scale of mu'
    std::optional<double> imaginaryOffset;  // log10 scale of mu''
    const CatalogueMaterial* donor;
};

// The donor's best frequency shift against the target's channels: at each shift each channel gets
// its own log offset (the amplitude scale) and the residual is the rms of what the offsets leave.
std::optional<DonorMatch> match_donor(const CatalogueMaterial& donor, const Channels& channels) {
    std::optional<DonorMatch> best;
    std::vector<double> residuals;
    std::vector<double> residualWeights;
    for (int shiftIndex = 0; shiftIndex < kNumberShifts; ++shiftIndex) {
        double shift = -kShiftRangeDecades + shiftIndex * kShiftStepDecades;
        residuals.clear();
        residualWeights.clear();
        std::optional<double> offsets[2];
        bool rejected = false;
        bool anyChannel = false;
        const Channel* channelList[2] = {&channels.magnitude, &channels.imaginary};
        for (size_t channelIndex = 0; channelIndex < 2 && !rejected; ++channelIndex) {
            const Channel& channel = *channelList[channelIndex];
            if (channel.logFrequencies.empty()) {
                continue;
            }
            bool isImaginary = channelIndex == 1;
            std::vector<double> differences;
            std::vector<double> weights;
            double totalWeight = 0;
            double matchedWeight = 0;
            for (size_t index = 0; index < channel.logFrequencies.size(); ++index) {
                totalWeight += channel.weights[index];
                auto donorValue = interpolate(isImaginary ? donor.imaginary : donor.real, channel.logFrequencies[index] - shift, isImaginary);
                if (donorValue && donorValue.value() > 0) {
                    differences.push_back(channel.logValues[index] - std::log10(donorValue.value()));
                    weights.push_back(channel.weights[index]);
                    matchedWeight += channel.weights[index];
                }
            }
            if (differences.size() < kMinimumMatchedPoints || matchedWeight < kMinimumMatchedFraction * totalWeight) {
                rejected = true;
                break;
            }
            double mean = 0;
            for (size_t index = 0; index < differences.size(); ++index) {
                mean += weights[index] * differences[index];
            }
            mean /= matchedWeight;
            offsets[channelIndex] = mean;
            for (size_t index = 0; index < differences.size(); ++index) {
                residuals.push_back(differences[index] - mean);
                residualWeights.push_back(weights[index]);
            }
            anyChannel = true;
        }
        if (rejected || !anyChannel) {
            continue;
        }
        double sumSquares = 0;
        double sumWeights = 0;
        for (size_t index = 0; index < residuals.size(); ++index) {
            sumSquares += residualWeights[index] * residuals[index] * residuals[index];
            sumWeights += residualWeights[index];
        }
        double residual = std::sqrt(sumSquares / sumWeights);
        if (!best || residual < best->residual) {
            best = DonorMatch{residual, shift, offsets[0], offsets[1], &donor};
        }
    }
    return best;
}

// Learned from the catalogue: the median log-log slope of the donors' outermost decade on that side
// (mu': of |mu'|; the sign of the estimate is kept). Disengaged when no donor has such a decade.
std::optional<std::pair<double, size_t>> learned_tail_slope(const std::vector<const CatalogueMaterial*>& donors,
                                                            ComplexPermeabilityPart part, ComplexPermeabilityGapSide side) {
    std::vector<double> slopes;
    for (const auto* donor : donors) {
        const auto& table = part == ComplexPermeabilityPart::IMAGINARY ? donor->imaginary : donor->real;
        const auto& x = table.logFrequencies;
        std::vector<double> xs, ys;
        bool positive = true;
        for (size_t index = 0; index < x.size(); ++index) {
            bool selected = side == ComplexPermeabilityGapSide::BELOW ? x[index] <= x.front() + 1.0 : x[index] >= x.back() - 1.0;
            if (!selected) {
                continue;
            }
            double value = std::fabs(table.values[index]);
            if (!(value > 0)) {
                positive = false;
                break;
            }
            xs.push_back(x[index]);
            ys.push_back(std::log10(value));
        }
        if (!positive || xs.size() < kTailMinimumPoints || !(xs.back() - xs.front() > std::log10(kTailMinimumSpan))) {
            continue;
        }
        slopes.push_back(linear_fit(xs, ys).first);
    }
    if (slopes.empty()) {
        return std::nullopt;
    }
    return std::make_pair(median(slopes), slopes.size());
}

// The median of the first kNeighbours donors (in match order) whose shifted table reaches the
// frequency; disengaged when fewer than kMinimumNeighbours do.
std::optional<double> neighbours_value(const std::vector<ComplexPermeabilityDonorFit>& donors, ComplexPermeabilityPart part, double logFrequency) {
    std::vector<double> values;
    bool logValues = part == ComplexPermeabilityPart::IMAGINARY;
    for (const auto& donor : donors) {
        auto value = interpolate(donor.table, logFrequency - donor.shiftDecades, logValues);
        if (value) {
            values.push_back(donor.gain * value.value());
            if (values.size() == kNeighbours) {
                break;
            }
        }
    }
    if (values.size() < kMinimumNeighbours) {
        return std::nullopt;
    }
    return median(values);
}

size_t donors_reaching(const std::vector<ComplexPermeabilityDonorFit>& donors, double logFrequency) {
    size_t count = 0;
    for (const auto& donor : donors) {
        const auto& x = donor.table.logFrequencies;
        double shifted = logFrequency - donor.shiftDecades;
        if (shifted >= x.front() - kLogTolerance && shifted <= x.back() + kLogTolerance) {
            ++count;
        }
    }
    return count;
}

// The outermost log-frequency, going away from the edge, up to which at least kMinimumNeighbours
// donors reach every point. Disengaged when they do not reach past the edge at all.
std::optional<double> neighbours_reach(const std::vector<ComplexPermeabilityDonorFit>& donors, double logEdge, ComplexPermeabilityGapSide side) {
    std::vector<double> boundaries;
    for (const auto& donor : donors) {
        for (double end : {donor.table.logFrequencies.front() + donor.shiftDecades, donor.table.logFrequencies.back() + donor.shiftDecades}) {
            if (beyond_edge(end, logEdge, side)) {
                boundaries.push_back(end);
            }
        }
    }
    std::sort(boundaries.begin(), boundaries.end());
    if (side == ComplexPermeabilityGapSide::BELOW) {
        std::reverse(boundaries.begin(), boundaries.end());
    }
    boundaries.erase(std::unique(boundaries.begin(), boundaries.end()), boundaries.end());
    std::optional<double> reach;
    double previous = logEdge;
    for (double boundary : boundaries) {
        if (donors_reaching(donors, (previous + boundary) / 2) < kMinimumNeighbours) {
            break;
        }
        reach = boundary;
        previous = boundary;
    }
    return reach;
}

// ABT #1870: mu' above its data from the material's OWN tabulated mu'' (never an estimated mu''),
// anchored at its last tabulated mu' point. Disengaged when its mu'' table does not reach past the mu'
// edge.
std::optional<ComplexPermeabilityKramersKronig> kramers_kronig_real(const CatalogueMaterial& target) {
    const auto& imaginary = target.imaginary;
    const auto& lf = imaginary.logFrequencies;
    const auto& lv = imaginary.logValues;
    double logRealEdge = target.real.logFrequencies.back();
    if (!(lf.back() > logRealEdge + kLogTolerance)) {
        return std::nullopt;
    }
    ComplexPermeabilityKramersKronig kk;
    kk.anchorFrequency = std::pow(10.0, logRealEdge);
    kk.anchorValue = target.real.values.back();
    kk.toFrequency = std::pow(10.0, lf.back());

    // The grid: numpy arange(lf[0] - pad, lf[-1] + pad, 1 / 200).
    const double step = 1.0 / kKramersKronigGridPerDecade;
    const double start = lf.front() - kKramersKronigPadDecades;
    const double stop = lf.back() + kKramersKronigPadDecades;
    const size_t numberPoints = static_cast<size_t>(std::ceil((stop - start) / step));
    std::vector<double> grid(numberPoints);
    std::vector<double> logGridValues(numberPoints);
    for (size_t index = 0; index < numberPoints; ++index) {
        grid[index] = start + static_cast<double>(index) * step;
        logGridValues[index] = interpolate_clamped(lf, lv, grid[index]);
    }
    // The pads: the power law of the table's outermost decade on each side (where it has 2 points).
    for (bool high : {false, true}) {
        std::vector<double> xs, ys;
        for (size_t index = 0; index < lf.size(); ++index) {
            bool selected = high ? lf[index] >= lf.back() - kKramersKronigEdgeFitDecades : lf[index] <= lf.front() + kKramersKronigEdgeFitDecades;
            if (selected) {
                xs.push_back(lf[index]);
                ys.push_back(lv[index]);
            }
        }
        if (xs.size() < 2) {
            continue;
        }
        double slope = linear_fit(xs, ys).first;
        if (high) {
            slope = std::min(slope, kKramersKronigMaximumHighSlope);
        }
        // The pad starts from the table's own end value. (kk.py keeps the fitted intercept: the pad then
        // jumps at the table's end wherever the slope was capped or the fit misses the end point, and the
        // transform has a logarithmic singularity there. Hold-out: identical MnZn cells, NiZn mu' P90 21 -> 7 %;
        // production: A13's mu' went from -1398 to -4861 within 3 % of frequency, R10KC's from 112 to -5904 within 6 %.)
        double intercept = high ? lv.back() - slope * lf.back() : lv.front() - slope * lf.front();
        for (size_t index = 0; index < numberPoints; ++index) {
            if (high ? grid[index] > lf.back() : grid[index] < lf.front()) {
                logGridValues[index] = slope * grid[index] + intercept;
            }
        }
    }
    kk.gridLogFrequencies.resize(numberPoints);
    kk.gridFrequencies.resize(numberPoints);
    kk.gridImaginary.resize(numberPoints);
    for (size_t index = 0; index < numberPoints; ++index) {
        kk.gridFrequencies[index] = std::pow(10.0, grid[index]);
        kk.gridLogFrequencies[index] = std::log(kk.gridFrequencies[index]);
        kk.gridImaginary[index] = std::pow(10.0, logGridValues[index]);
    }
    kk.transformAtAnchor = kk.transform(kk.anchorFrequency);
    std::ostringstream method;
    method << "mu' ESTIMATED by Kramers-Kronig from " << target.name << "'s own tabulated mu'' (" << format_frequency(std::pow(10.0, lf.front()))
           << " to " << format_frequency(kk.toFrequency) << "), anchored at its last tabulated mu', " << kk.anchorValue << " at "
           << format_frequency(kk.anchorFrequency);
    kk.method = method.str();
    return kk;
}

std::string class_name(const CoreMaterial& material) {
    return std::string(magic_enum::enum_name(material.get_material())) + " " +
           std::string(magic_enum::enum_name(material.get_material_composition().value()));
}

ComplexPermeabilityGapEstimate held_flat_gap(const CatalogueMaterial& target) {
    ComplexPermeabilityGapEstimate gap;
    gap.materialName = target.name;
    gap.part = ComplexPermeabilityPart::REAL;
    gap.side = ComplexPermeabilityGapSide::BELOW;
    gap.edgeFrequency = std::pow(10.0, target.real.logFrequencies.front());
    gap.source = ComplexPermeabilitySource::HELD_FLAT;
    gap.heldValue = target.real.values.front();
    std::ostringstream method;
    method << "mu' held at its lowest tabulated value, " << gap.heldValue << " at " << format_frequency(gap.edgeFrequency)
           << " (the catalogue's mu' is flat below its relaxation)";
    gap.neighboursMethod = method.str();
    return gap;
}

ComplexPermeabilityGapEstimate neighbours_gap(const CatalogueMaterial& target, const std::string& materialClass,
                                              const std::vector<const CatalogueMaterial*>& donors,
                                              ComplexPermeabilityPart part, ComplexPermeabilityGapSide side) {
    ComplexPermeabilityGapEstimate gap;
    gap.materialName = target.name;
    gap.part = part;
    gap.side = side;
    gap.source = ComplexPermeabilitySource::ESTIMATED_NEIGHBOURS;
    const auto& own = part == ComplexPermeabilityPart::REAL ? target.real : target.imaginary;
    const auto& other = part == ComplexPermeabilityPart::REAL ? target.imaginary : target.real;
    double logEdge = side == ComplexPermeabilityGapSide::BELOW ? own.logFrequencies.front() : own.logFrequencies.back();
    gap.edgeFrequency = std::pow(10.0, logEdge);
    // The other part is known inside the gap when its table reaches past this edge.
    bool otherKnownBeyond = side == ComplexPermeabilityGapSide::BELOW ? other.logFrequencies.front() < logEdge
                                                                       : other.logFrequencies.back() > logEdge;
    bool realMissing = part == ComplexPermeabilityPart::REAL || !otherKnownBeyond;
    bool imaginaryMissing = part == ComplexPermeabilityPart::IMAGINARY || !otherKnownBeyond;
    auto channels = match_channels(target, logEdge, side, realMissing, imaginaryMissing);
    std::string where = part_symbol(part) + (side == ComplexPermeabilityGapSide::BELOW ? " below " : " above ") + format_frequency(gap.edgeFrequency);
    if (channels.magnitude.logFrequencies.size() + channels.imaginary.logFrequencies.size() < kMinimumMatchedPoints) {
        gap.refusal = where + ": fewer than " + std::to_string(kMinimumMatchedPoints) + " of its own points to match the donors against";
        return gap;
    }

    std::vector<DonorMatch> matches;
    for (const auto* donor : donors) {
        if (auto match = match_donor(*donor, channels)) {
            matches.push_back(match.value());
        }
    }
    std::stable_sort(matches.begin(), matches.end(), [](const DonorMatch& a, const DonorMatch& b) { return a.residual < b.residual; });
    gap.numberDonorsMatched = matches.size();
    if (matches.size() < kMinimumNeighbours) {
        gap.refusal = where + ": only " + std::to_string(matches.size()) + " " + materialClass + " donors match its own data (at least " +
                      std::to_string(kMinimumNeighbours) + " needed)";
        return gap;
    }
    for (size_t index = 0; index < std::min(kDonorsConsidered, matches.size()); ++index) {
        const auto& match = matches[index];
        ComplexPermeabilityDonorFit fit;
        fit.name = match.donor->name;
        fit.manufacturer = match.donor->manufacturer;
        fit.shiftDecades = match.shiftDecades;
        fit.residual = match.residual;
        // mu' takes the mu' (magnitude) channel's scale, mu'' the mu'' channel's; each falls back on
        // the other channel's when its own was not matched (no points of it next to the edge).
        double offset = part == ComplexPermeabilityPart::REAL ? match.magnitudeOffset.value_or(match.imaginaryOffset.value_or(0))
                                                              : match.imaginaryOffset.value_or(match.magnitudeOffset.value_or(0));
        fit.gain = std::pow(10.0, offset);
        fit.table = part == ComplexPermeabilityPart::REAL ? match.donor->real : match.donor->imaginary;
        gap.donors.push_back(std::move(fit));
    }

    auto reach = neighbours_reach(gap.donors, logEdge, side);
    if (!reach) {
        gap.refusal = where + ": fewer than " + std::to_string(kMinimumNeighbours) + " of its best-matching " + materialClass +
                      " donors have data past the edge";
        return gap;
    }
    // The donors give the shape; the estimate meets the material's own last tabulated point exactly
    // (Alf, 2026-10-14). mu'': every donor's gain is multiplied by own edge value / donors' median
    // there. mu': the difference is added, fading over kEdgeOffsetFadeDecades (edge_offset_at), so a
    // donors' median on the other side of zero is joined, not flipped.
    {
        auto edgeEstimate = neighbours_value(gap.donors, part, gap.inside_gap(logEdge));
        if (!edgeEstimate) {
            gap.refusal = where + ": fewer than " + std::to_string(kMinimumNeighbours) + " of its best-matching " + materialClass +
                          " donors reach the edge, so the estimate cannot be joined to its own value there";
            return gap;
        }
        if (edgeEstimate.value() == 0) {
            gap.refusal = where + ": the donors' median at the edge is 0 (their tables are 0 there), so they give no trend to join "
                          "its own value to";
            return gap;
        }
        double ownEdge = side == ComplexPermeabilityGapSide::BELOW ? own.values.front() : own.values.back();
        if (part == ComplexPermeabilityPart::REAL) {
            gap.edgeOffset = ownEdge - edgeEstimate.value();
            gap.edgeSignMismatch = ownEdge * edgeEstimate.value() < 0;
        }
        else {
            gap.edgeScale = ownEdge / edgeEstimate.value();
            for (auto& fit : gap.donors) {
                fit.gain *= gap.edgeScale;
            }
        }
    }
    gap.reachFrequency = std::pow(10.0, reach.value());
    auto reachValue = neighbours_value(gap.donors, part, reach.value());
    if (!reachValue) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
            "Complex permeability estimate of " + target.name + ", " + where + ": the donors' reach has no estimate");
    }
    // mu': the donors' value, without the edge offset; evaluate() adds the offset in the tail as well,
    // so it is gone a decade from the edge wherever that falls.
    gap.reachValue = reachValue.value();
    auto slope = learned_tail_slope(donors, part, side);
    if (!slope) {
        gap.refusal = where + ": no " + materialClass + " donor has an outermost decade to learn the tail slope from";
        return gap;
    }
    gap.tailSlope = slope->first;
    gap.tailSlopeDonors = slope->second;

    std::ostringstream neighbours;
    neighbours << part_symbol(part) << " ESTIMATED from the " << std::min(kNeighbours, gap.donors.size()) << " best shape-matched "
               << materialClass << " donors (of " << gap.numberDonorsMatched << " matched):";
    for (size_t index = 0; index < std::min(kNeighbours, gap.donors.size()); ++index) {
        neighbours << (index == 0 ? " " : ", ") << gap.donors[index].name << " " << std::showpos << std::fixed << std::setprecision(3)
                   << gap.donors[index].shiftDecades << std::noshowpos << " dec";
    }
    neighbours << "; per-point median, each donor shifted in frequency and scaled to " << target.name << "'s own data ("
               << (side == ComplexPermeabilityGapSide::ABOVE ? "weighted down with distance from " : "within a decade of ")
               << format_frequency(gap.edgeFrequency) << "), then ";
    if (part == ComplexPermeabilityPart::REAL) {
        neighbours << "offset by " << std::showpos << std::defaultfloat << std::setprecision(4) << gap.edgeOffset << std::noshowpos
                   << " at " << format_frequency(gap.edgeFrequency) << ", fading to 0 a decade from it, to meet its own value there";
    }
    else {
        neighbours << "scaled by " << std::defaultfloat << std::setprecision(4) << gap.edgeScale << " to meet its own value at "
                   << format_frequency(gap.edgeFrequency);
    }
    neighbours << "; reaches " << format_frequency(gap.reachFrequency);
    gap.neighboursMethod = neighbours.str();
    std::ostringstream tail;
    tail << part_symbol(part) << " ESTIMATED beyond the donors' reach (" << format_frequency(gap.reachFrequency) << ") along the "
         << materialClass << " tail slope " << std::showpos << std::fixed << std::setprecision(3) << gap.tailSlope << std::noshowpos
         << " per decade (median over " << gap.tailSlopeDonors << " donors)";
    gap.tailMethod = tail.str();
    return gap;
}

}  // namespace

double ComplexPermeabilityKramersKronig::transform(double frequency) const {
    // Port of kk.py kk_real: (vg fg^2 - m0 f0^2) / (fg^2 - f0^2) integrated by the trapezoid rule in
    // ln f over the grid points away from the singularity, plus the analytic principal value of
    // m0 f0^2 / (fg^2 - f0^2) over the grid, (m0 / 2) [ln|1 - f0^2 / fg^2|] between the grid's ends.
    const double f0Squared = frequency * frequency;
    const double m0 = interpolate_clamped(gridLogFrequencies, gridImaginary, std::log(frequency));
    double integral = 0;
    bool havePrevious = false;
    double previousU = 0, previousValue = 0;
    for (size_t index = 0; index < gridFrequencies.size(); ++index) {
        double fgSquared = gridFrequencies[index] * gridFrequencies[index];
        double denominator = fgSquared - f0Squared;
        if (!(std::fabs(denominator) > kKramersKronigSingularTolerance * f0Squared)) {
            continue;
        }
        double value = (gridImaginary[index] * fgSquared - m0 * f0Squared) / denominator;
        double u = gridLogFrequencies[index];
        if (havePrevious) {
            integral += (u - previousU) * (value + previousValue) / 2;
        }
        previousU = u;
        previousValue = value;
        havePrevious = true;
    }
    auto principal = [&](double gridFrequency) {
        return std::log(std::fabs(1 - f0Squared / (gridFrequency * gridFrequency)));
    };
    double singular = m0 / 2 * (principal(gridFrequencies.back()) - principal(gridFrequencies.front()));
    return 1 + 2 / std::numbers::pi * (integral + singular);
}

double ComplexPermeabilityKramersKronig::real_part(double frequency) const {
    return anchorValue + transform(frequency) - transformAtAnchor;
}

std::pair<double, ComplexPermeabilitySource> ComplexPermeabilityGapEstimate::evaluate(double frequency) const {
    // ABT #1870 (NiZn only): mu' where the own mu'' is tabulated comes from it, whatever the neighbours
    // say (or refuse); over the last 0.25 decade of that table it moves linearly in log f into the
    // neighbours, which alone continue beyond it. Without KK this is the neighbours exactly (ABT #1861).
    if (kramersKronig) {
        if (!std::isfinite(frequency) || frequency <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Complex permeability of " + materialName + " asked at " +
                                        std::to_string(frequency) + " Hz: a frequency must be a finite positive number");
        }
        if (frequency <= kramersKronig->toFrequency) {
            double kramersKronigValue = kramersKronig->real_part(frequency);
            if (refusal.empty()) {
                double start = kramers_kronig_blend_start();
                if (frequency > start) {
                    double weight = std::log10(frequency / start) / std::log10(kramersKronig->toFrequency / start);
                    return {(1 - weight) * kramersKronigValue + weight * evaluate_neighbours(frequency).first,
                            ComplexPermeabilitySource::ESTIMATED_KRAMERS_KRONIG};
                }
            }
            return {kramersKronigValue, ComplexPermeabilitySource::ESTIMATED_KRAMERS_KRONIG};
        }
    }
    return evaluate_neighbours(frequency);
}

double ComplexPermeabilityGapEstimate::kramers_kronig_blend_start() const {
    if (!kramersKronig) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Complex permeability estimate of " + materialName +
                                   ": a Kramers-Kronig blend was asked for a gap without a Kramers-Kronig range");
    }
    return std::max(edgeFrequency, kramersKronig->toFrequency / std::pow(10.0, kKramersKronigBlendDecades));
}

std::pair<double, ComplexPermeabilitySource> ComplexPermeabilityGapEstimate::evaluate_neighbours(double frequency) const {
    if (!refusal.empty()) {
        throw ComplexPermeabilityEstimatorRefusal(materialName, refusal);
    }
    if (!std::isfinite(frequency) || frequency <= 0) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Complex permeability of " + materialName + " asked at " +
                                    std::to_string(frequency) + " Hz: a frequency must be a finite positive number");
    }
    if (source == ComplexPermeabilitySource::HELD_FLAT) {
        return {heldValue, ComplexPermeabilitySource::HELD_FLAT};
    }
    double logFrequency = std::log10(frequency);
    double logReach = std::log10(reachFrequency);
    bool beyondReach = side == ComplexPermeabilityGapSide::BELOW ? logFrequency < logReach : logFrequency > logReach;
    if (beyondReach) {
        return {reachValue * std::pow(frequency / reachFrequency, tailSlope) + edge_offset_at(frequency), ComplexPermeabilitySource::ESTIMATED_TAIL};
    }
    auto value = neighbours_value(donors, part, inside_gap(logFrequency));
    if (!value) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
            "Complex permeability estimate of " + materialName + ": fewer than 3 donors reach " + std::to_string(frequency) +
            " Hz inside their reach (" + std::to_string(reachFrequency) + " Hz)");
    }
    return {value.value() + edge_offset_at(frequency), ComplexPermeabilitySource::ESTIMATED_NEIGHBOURS};
}

std::optional<double> ComplexPermeabilityGapEstimate::neighbours_median(double frequency) const {
    if (source != ComplexPermeabilitySource::ESTIMATED_NEIGHBOURS || !refusal.empty()) {
        return std::nullopt;
    }
    return neighbours_value(donors, part, inside_gap(std::log10(frequency)));
}

double ComplexPermeabilityGapEstimate::inside_gap(double logFrequency) const {
    double logEdge = std::log10(edgeFrequency);
    return side == ComplexPermeabilityGapSide::ABOVE ? std::max(logFrequency, logEdge + kEdgeJoinStepDecades)
                                                     : std::min(logFrequency, logEdge - kEdgeJoinStepDecades);
}

double ComplexPermeabilityGapEstimate::edge_offset_at(double frequency) const {
    if (part != ComplexPermeabilityPart::REAL || source != ComplexPermeabilitySource::ESTIMATED_NEIGHBOURS) {
        return 0;
    }
    double distance = std::fabs(std::log10(frequency / edgeFrequency));
    return edgeOffset * std::max(0.0, 1 - distance / kEdgeOffsetFadeDecades);
}

const std::string& ComplexPermeabilityGapEstimate::method(ComplexPermeabilitySource valueSource) const {
    if (valueSource == ComplexPermeabilitySource::ESTIMATED_KRAMERS_KRONIG) {
        if (!kramersKronig) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Complex permeability estimate of " + materialName +
                                       ": a Kramers-Kronig value was asked for a gap without a Kramers-Kronig range");
        }
        return kramersKronig->method;
    }
    return valueSource == ComplexPermeabilitySource::ESTIMATED_TAIL ? tailMethod : neighboursMethod;
}

const ComplexPermeabilityGapEstimate& ComplexPermeabilityEstimate::gap(ComplexPermeabilityPart part, ComplexPermeabilityGapSide side) const {
    for (const auto& candidate : gaps) {
        if (candidate.part == part && candidate.side == side) {
            return candidate;
        }
    }
    throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Complex permeability estimate of " + materialName + " has no " +
                               part_symbol(part) + " gap on that side");
}

std::string ComplexPermeabilityEstimate::summary() const {
    std::ostringstream text;
    text << "Core material '" << materialName << "': complex permeability is tabulated for mu' from "
         << format_frequency(realSpan.first) << " to " << format_frequency(realSpan.second) << " and for mu'' from "
         << format_frequency(imaginarySpan.first) << " to " << format_frequency(imaginarySpan.second)
         << "; outside that it is ESTIMATED from " << numberClassDonors << " " << materialClass << " donors";
    if (!snoekScreened.empty()) {
        text << " (" << snoekScreened.size() << " screened out: mu_i * f(mu'' peak) more than 2 decades off the class median)";
    }
    text << ".";
    for (const auto& gap : gaps) {
        text << " " << (gap.side == ComplexPermeabilityGapSide::BELOW ? "Below " : "Above ") << format_frequency(gap.edgeFrequency) << ": ";
        if (gap.kramersKronig) {
            text << gap.kramersKronig->method << " up to " << format_frequency(gap.kramersKronig->toFrequency) << "; beyond it: ";
        }
        if (!gap.refusal.empty()) {
            text << "NOT estimated (" << gap.refusal << ").";
        }
        else if (gap.source == ComplexPermeabilitySource::HELD_FLAT) {
            text << gap.neighboursMethod << ".";
        }
        else {
            text << gap.neighboursMethod << "; " << gap.tailMethod << ".";
        }
    }
    return text.str();
}

ComplexPermeabilityEstimate estimate_complex_permeability(const CoreMaterial& coreMaterial) {
    const std::string& materialName = coreMaterial.get_name();
    const auto composition = coreMaterial.get_material_composition();
    if (!composition) {
        throw ComplexPermeabilityEstimatorRefusal(materialName, "it has no materialComposition, and the estimator only covers the "
                                                  "FERRITE MN_ZN and FERRITE NI_ZN classes");
    }
    if (coreMaterial.get_material() != MaterialType::FERRITE ||
        (composition.value() != MaterialComposition::MN_ZN && composition.value() != MaterialComposition::NI_ZN)) {
        throw ComplexPermeabilityEstimatorRefusal(materialName, "its class " + class_name(coreMaterial) + " is not covered: the estimator "
                                                  "only covers FERRITE MN_ZN and FERRITE NI_ZN (powders, MgZn and other classes are too "
                                                  "few and too heterogeneous in MAS to estimate from, proposal hold-out)");
    }
    if (!coreMaterial.get_permeability().get_complex()) {
        throw ComplexPermeabilityEstimatorRefusal(materialName, "it has no tabulated complex permeability in MAS (a table derived from "
                                                  "its initial permeability is not estimated beyond its range)");
    }
    auto target = catalogue_material(coreMaterial);
    if (!target) {
        throw ComplexPermeabilityEstimatorRefusal(materialName, "its tabulated mu' or mu'' has fewer than 2 usable points");
    }
    const std::string materialClass = class_name(coreMaterial);

    if (coreMaterialDatabase.empty()) {
        load_core_materials();
    }
    std::vector<CatalogueMaterial> members;
    for (const auto& [name, material] : coreMaterialDatabase) {
        if (material.get_material() != coreMaterial.get_material() || material.get_material_composition() != composition) {
            continue;
        }
        if (auto entry = catalogue_material(material)) {
            members.push_back(std::move(entry.value()));
        }
    }

    // Snoek screen, learned from the class: K = median of mu_i * f(mu'' peak).
    std::vector<double> snoekProducts;
    for (const auto& member : members) {
        if (member.lossPeakFrequency) {
            snoekProducts.push_back(member.initialPermeability * member.lossPeakFrequency.value());
        }
    }
    ComplexPermeabilityEstimate estimate;
    estimate.materialName = materialName;
    estimate.materialClass = materialClass;
    std::optional<double> snoekProduct;
    if (snoekProducts.size() >= kSnoekScreenMinimumMembers) {
        snoekProduct = median(snoekProducts);
    }
    std::vector<const CatalogueMaterial*> donors;
    for (const auto& member : members) {
        if (snoekProduct && member.lossPeakFrequency &&
            std::fabs(std::log10(member.initialPermeability * member.lossPeakFrequency.value() / snoekProduct.value())) > kSnoekScreenDecades) {
            estimate.snoekScreened.push_back(member.name);
            continue;
        }
        if (member.name == target->name || is_near_duplicate(target.value(), member)) {
            continue;
        }
        donors.push_back(&member);
    }
    estimate.numberClassDonors = donors.size();
    if (donors.size() < kMinimumNeighbours) {
        throw ComplexPermeabilityEstimatorRefusal(materialName, "its class " + materialClass + " has only " + std::to_string(donors.size()) +
                                                  " donors with tabulated complex permeability in MAS (at least " +
                                                  std::to_string(kMinimumNeighbours) + " needed)");
    }

    auto rawSpan = [&](const MAS::Permeability& part) {
        double minimum = std::numeric_limits<double>::max();
        double maximum = std::numeric_limits<double>::lowest();
        for (const auto& point : tabulated_points(part)) {
            minimum = std::min(minimum, point.get_frequency().value());
            maximum = std::max(maximum, point.get_frequency().value());
        }
        return std::make_pair(minimum, maximum);
    };
    estimate.realSpan = rawSpan(coreMaterial.get_permeability().get_complex()->get_real());
    estimate.imaginarySpan = rawSpan(coreMaterial.get_permeability().get_complex()->get_imaginary());

    estimate.gaps.push_back(held_flat_gap(target.value()));
    estimate.gaps.push_back(neighbours_gap(target.value(), materialClass, donors, ComplexPermeabilityPart::REAL, ComplexPermeabilityGapSide::ABOVE));
    // ABT #1870: a NiZn ferrite's mu' above its data comes from its own tabulated mu'' where that reaches
    // further (anchored KK); the neighbours take over over the last 0.25 decade of the mu'' table and
    // continue alone beyond it. MnZn keeps the neighbours (kKramersKronigComposition).
    if (composition.value() == kKramersKronigComposition) {
        if (auto kramersKronig = kramers_kronig_real(target.value())) {
            auto& realAbove = estimate.gaps.back();
            realAbove.kramersKronig = std::move(kramersKronig);
            std::ostringstream join;
            if (realAbove.refusal.empty()) {
                join << "; from " << format_frequency(realAbove.kramers_kronig_blend_start()) << " to "
                     << format_frequency(realAbove.kramersKronig->toFrequency)
                     << " (the last 0.25 decade of its own mu'') moved linearly in log f into the neighbours' mu', which continue alone beyond it";
            }
            else {
                join << "; beyond its own mu'' it is NOT estimated (" << realAbove.refusal << ")";
            }
            realAbove.kramersKronig->method += join.str();
        }
    }
    estimate.gaps.push_back(neighbours_gap(target.value(), materialClass, donors, ComplexPermeabilityPart::IMAGINARY, ComplexPermeabilityGapSide::BELOW));
    estimate.gaps.push_back(neighbours_gap(target.value(), materialClass, donors, ComplexPermeabilityPart::IMAGINARY, ComplexPermeabilityGapSide::ABOVE));
    return estimate;
}

std::vector<ComplexPermeabilityEstimatedRange> get_estimated_ranges(const ComplexPermeabilityEstimate& estimate) {
    std::vector<ComplexPermeabilityEstimatedRange> ranges;
    const double infinity = std::numeric_limits<double>::infinity();
    for (const auto& gap : estimate.gaps) {
        const auto& span = gap.part == ComplexPermeabilityPart::REAL ? estimate.realSpan : estimate.imaginarySpan;
        double edge = gap.side == ComplexPermeabilityGapSide::BELOW ? span.first : span.second;
        if (gap.kramersKronig) {
            // ABT #1870: mu' from its own tabulated mu'' first; the neighbours' ranges start where it ends.
            ranges.push_back({edge, gap.kramersKronig->toFrequency, gap.part, ComplexPermeabilitySource::ESTIMATED_KRAMERS_KRONIG,
                              estimate.materialClass, {}, gap.kramersKronig->method});
            edge = gap.kramersKronig->toFrequency;
        }
        if (!gap.refusal.empty()) {
            continue;
        }
        std::vector<std::string> donorNames;
        for (size_t index = 0; index < std::min(kNeighbours, gap.donors.size()); ++index) {
            donorNames.push_back(gap.donors[index].name);
        }
        if (gap.source == ComplexPermeabilitySource::HELD_FLAT) {
            ranges.push_back({0, edge, gap.part, gap.source, estimate.materialClass, {}, gap.neighboursMethod});
            continue;
        }
        if (gap.side == ComplexPermeabilityGapSide::BELOW) {
            ranges.push_back({gap.reachFrequency, edge, gap.part, ComplexPermeabilitySource::ESTIMATED_NEIGHBOURS, estimate.materialClass, donorNames, gap.neighboursMethod});
            ranges.push_back({0, gap.reachFrequency, gap.part, ComplexPermeabilitySource::ESTIMATED_TAIL, estimate.materialClass, donorNames, gap.tailMethod});
        }
        else {
            // With KK the neighbours' ranges start at the end of the own mu'' table, which may lie past their reach.
            if (!gap.kramersKronig || gap.reachFrequency > edge) {
                ranges.push_back({edge, gap.reachFrequency, gap.part, ComplexPermeabilitySource::ESTIMATED_NEIGHBOURS, estimate.materialClass, donorNames, gap.neighboursMethod});
            }
            ranges.push_back({gap.kramersKronig ? std::max(edge, gap.reachFrequency) : gap.reachFrequency, infinity, gap.part, ComplexPermeabilitySource::ESTIMATED_TAIL, estimate.materialClass, donorNames, gap.tailMethod});
        }
    }
    return ranges;
}

}  // namespace OpenMagnetics
