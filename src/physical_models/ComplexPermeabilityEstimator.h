#pragma once
#include "support/Exceptions.h"
#include <MAS.hpp>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

using namespace MAS;

namespace OpenMagnetics {

/**
 * THE APPROVED ESTIMATE of a core material's complex permeability outside its tabulated range
 * (ABT #1861, Alf's decision 2026-10-14). It is the second sanctioned estimate in MKF, after the
 * heat-conductivity estimator (resolve_core_material_heat_conductivity, ABT #1752/#1753), and it
 * follows the same rules: every number comes from the MAS catalogue at call time, nothing physical
 * is hardcoded, every estimated value is labelled with where it came from, and a material is
 * warned about once.
 *
 * Alf, 2026-10-14: complex permeability always extrapolates (and only complex permeability:
 * Settings::allowMaterialDataExtrapolation stays opt-in for every other material property); the
 * extrapolation estimates a material's missing range from the other materials' data; the impedance
 * sweeps return the extrapolated ranges next to the curve.
 *
 * The method (proposal .cperm-estimator/PROPOSAL.md, validated by hold-out on the catalogue):
 * - Shape-matched neighbours. Every other material of the same class (FERRITE MN_ZN, or FERRITE
 *   NI_ZN; same manufacturer allowed) is slid along log-frequency (the Snoek degree of freedom) and
 *   its mu' and mu'' are scaled, to fit the target's OWN data on the known side of the gap, plus the
 *   other part inside the gap, where that part is tabulated. Above the data every own point counts,
 *   weighted down with distance from the edge (exp(-d / 1 decade)); below it, the own points within
 *   one decade of the edge count, uniformly (Alf, 2026-10-14: each side uses the variant that did best
 *   on its own hold-out). The estimate at a frequency is the median of the 5 best-matching donors
 *   whose shifted data reach it (at least 3 required), joined to the material's own last tabulated
 *   value at the edge (the donors give the shape, the material the level): mu'' is scaled to meet
 *   it; mu', which crosses zero above resonance, takes the gap at the edge as an additive offset that
 *   fades linearly in log-frequency to zero one decade from the edge, so the curve is continuous
 *   there and follows the donors' trend beyond (a scale flipped the sign of mu' where the donors'
 *   median and the own value lie on opposite sides of zero).
 * - Beyond the donors' reach, the estimate continues from its outermost point along the class's tail
 *   slope, learned from the catalogue (median log-log slope of the donors' outermost decade).
 * - Below its data, mu' is held at its lowest tabulated value (the catalogue shows mu' flat there).
 * - Donors whose mu_i * f(mu'' peak) is more than 2 decades off their class median are screened out
 *   (19 DMEGC records and ACME N5 sit ~3 decades low: frequencies that look like kHz stored as Hz).
 * Hold-out errors of |Z| (toroid) with data removed above 1 MHz, MnZn: median 17 %, P90 66 % (C++ port,
 * with the edge join and the weights above; the proposal's Python: 19 % / 68 %), against
 * 182 % / 1696 % for the best simple baseline.
 *
 * Anchored Kramers-Kronig, FERRITE NI_ZN only (ABT #1870, Alf 2026-10-14): where a NiZn ferrite's mu' is
 * missing above its data but its OWN mu'' is tabulated further, mu' there is the singly-subtracted
 * Kramers-Kronig transform of that tabulated mu'', anchored at the last tabulated mu' point:
 *     mu'(f) = mu'(f_a) + KK[mu''](f) - KK[mu''](f_a),
 *     KK[mu''](f) = 1 + (2/pi) PV int mu''(f') f'^2 / (f'^2 - f^2) d ln f'
 * (port of .cperm-estimator/physical/kk.py and knn_kk.py: mu'' on a 200-per-decade log grid, padded one
 * decade each side with its edge decade's power-law slope, the high slope capped at -0.5, starting from
 * the table's end value (kk.py starts from the fitted line, which jumps at the end); the singularity
 * subtracted analytically). Estimated mu'' is NEVER fed into it (the research shows it makes |Z| worse):
 * the pad is the own table's power law, not the neighbours. Over the last 0.25 decade of the own mu''
 * table KK moves linearly in log f into the neighbours' mu', which continue alone beyond it.
 * Hold-out (TDK pattern, mu' hidden above the mu'' peak; C++ port): NiZn mu' median / P90 error 1 / 7 %
 * against 8 / 55 % for the neighbours. MnZn keeps the neighbours: there the anchor offset
 * mu'(f_a) - KK(f_a) is mu'(infinity) - 1, and on TDK MnZn records whose mu' and mu'' curves were
 * digitised slightly apart in frequency it sets implausible high-frequency levels (N51: +465 at 1 GHz)
 * that no check on the overlap can see. The datasheet lossFactor (tan d / mu_i) was tried as an anchor
 * for mu'' below the data and lost to the neighbours on the hold-out, so it is not used.
 *
 * It refuses (ComplexPermeabilityEstimatorRefusal) powders, MgZn, a proprietary or missing
 * composition, materials without tabulated complex data, a class with fewer than 3 donors, and a
 * target whose data next to the gap cannot be matched. The advisers and filters never use it: they
 * stay on the tabulated span (ComplexPermeability::get_tabulated_complex_permeability).
 */
enum class ComplexPermeabilitySource {
    TABULATED,             // the material's own MAS table
    ESTIMATED_NEIGHBOURS,  // median of the shape-matched class donors
    ESTIMATED_TAIL,        // beyond the donors' reach: the class tail slope
    HELD_FLAT,             // mu' below its data: its lowest tabulated value
    ESTIMATED_KRAMERS_KRONIG  // NiZn mu' above its data, from its own tabulated mu'' (ABT #1870)
};

enum class ComplexPermeabilityPart { REAL, IMAGINARY };

enum class ComplexPermeabilityGapSide { BELOW, ABOVE };

/**
 * @brief The estimator cannot estimate this material's missing range; reason says why.
 * ComplexPermeability turns it into ComplexPermeabilityEstimateUnavailableException, which names the
 * frequency that was asked for.
 */
class ComplexPermeabilityEstimatorRefusal : public MaterialException {
    std::string _reason;
public:
    ComplexPermeabilityEstimatorRefusal(const std::string& materialName, const std::string& reason)
        : MaterialException(ErrorCode::MATERIAL_DATA_MISSING,
                            "Complex permeability of " + materialName + " cannot be estimated: " + reason, materialName),
          _reason(reason) {}
    const std::string& reason() const noexcept { return _reason; }
};

// A cleaned complex-permeability part: sorted by frequency, exact duplicate frequencies averaged,
// mu'' keeping only its positive points.
struct ComplexPermeabilityTable {
    std::vector<double> logFrequencies;  // log10(Hz)
    std::vector<double> values;
    std::vector<double> logValues;       // log10(values), filled for mu'' only
};

// One donor kept for a gap: where it was shifted and how it was scaled to the target.
struct ComplexPermeabilityDonorFit {
    std::string name;
    std::string manufacturer;
    double shiftDecades;    // the donor is read at f / 10^shift
    double residual;        // rms of the log10 residual of the match
    double gain;            // scale applied to the donor's estimated part
    ComplexPermeabilityTable table;  // the donor's table of the estimated part
};

/**
 * @brief mu' above its tabulated data from the material's own tabulated mu'' by the anchored
 * Kramers-Kronig transform (ABT #1870). Built from the own mu'' table only.
 */
struct ComplexPermeabilityKramersKronig {
    double anchorFrequency;   // the last tabulated mu' point, Hz
    double anchorValue;       // mu' there
    double toFrequency;       // the last tabulated (positive) mu'' point, Hz: the transform is used up to here
    std::vector<double> gridLogFrequencies;  // ln f of the grid (natural log, as kk.py)
    std::vector<double> gridFrequencies;     // Hz
    std::vector<double> gridImaginary;       // the own mu'' on the grid (with the power-law pads)
    double transformAtAnchor = 0;            // KK[mu''](anchorFrequency)
    std::string method;

    // KK[mu''](f) = 1 + (2/pi) PV int mu''(f') f'^2 / (f'^2 - f^2) d ln f', singularity subtracted.
    double transform(double frequency) const;
    // mu'(f) = anchorValue + transform(f) - transformAtAnchor.
    double real_part(double frequency) const;
};

/**
 * @brief The estimate of one part (mu' or mu'') on one side of its tabulated data.
 */
struct ComplexPermeabilityGapEstimate {
    std::string materialName;
    ComplexPermeabilityPart part;
    ComplexPermeabilityGapSide side;
    double edgeFrequency;  // the part's last tabulated frequency on that side, Hz
    // HELD_FLAT (mu' below its data) or ESTIMATED_NEIGHBOURS (then ESTIMATED_TAIL past the reach).
    ComplexPermeabilitySource source;
    double heldValue = 0;                        // HELD_FLAT only
    std::vector<ComplexPermeabilityDonorFit> donors;  // ranked, the best first (ESTIMATED_NEIGHBOURS)
    size_t numberDonorsMatched = 0;              // class donors that could be matched at all
    double reachFrequency = 0;                   // outermost frequency at least 3 donors reach, Hz
    double reachValue = 0;                       // the donors' value there (mu': before the edge offset)
    double tailSlope = 0;                        // learned: d log|value| / d log f past the reach
    size_t tailSlopeDonors = 0;                  // donors the tail slope is the median of
    // The join to the material's own last tabulated value at the edge (Alf, 2026-10-14):
    // - mu'': every donor's gain carries edgeScale = own edge value / the donors' median there;
    // - mu': edgeOffset = own edge value - the donors' median there, added at the edge and fading
    //   linearly in log-frequency to 0 one decade from it (edge_offset_at). edgeScale stays 1.
    double edgeScale = 1;
    double edgeOffset = 0;
    // mu' only: the donors' median at the edge lies on the other side of zero from the own value. The
    // additive offset joins the two without flipping the donors' trend.
    bool edgeSignMismatch = false;
    // mu' above its data only, FERRITE NI_ZN only (ABT #1870): where the own mu'' is tabulated beyond the mu'
    // edge, mu' up to kramersKronig->toFrequency is the anchored KK transform (ESTIMATED_KRAMERS_KRONIG),
    // even when the neighbours are refused. Over the last 0.25 decade of the own mu'' table it moves
    // linearly in log f from KK to the neighbours (when they answer); beyond the table, the neighbours
    // alone, exactly as without KK. Never engaged for MnZn.
    std::optional<ComplexPermeabilityKramersKronig> kramersKronig;
    std::string neighboursMethod;                // human-readable, as CoreHeatConductivityEstimate::method
    std::string tailMethod;
    // Non-empty when this part and side cannot be estimated (e.g. no donor's data reach past the
    // edge): evaluate() throws ComplexPermeabilityEstimatorRefusal with it. The other gaps stand.
    std::string refusal;

    // The estimated value at a frequency beyond the edge, and where it came from.
    std::pair<double, ComplexPermeabilitySource> evaluate(double frequency) const;
    // The neighbours' (or their tail's) value, with the mu' edge offset (ABT #1861), without Kramers-Kronig.
    std::pair<double, ComplexPermeabilitySource> evaluate_neighbours(double frequency) const;
    // Where the blend from Kramers-Kronig into the neighbours starts, Hz (only with kramersKronig and
    // the neighbours answering): 0.25 decade below the end of the own mu'' table, not below the edge.
    double kramers_kronig_blend_start() const;
    // mu': the part of edgeOffset still added at this frequency (edgeOffset at the edge, 0 from one
    // decade beyond it). 0 for mu'' and for HELD_FLAT.
    double edge_offset_at(double frequency) const;
    // The donors' per-point median at a frequency inside their reach, before the mu' edge offset (for
    // mu'' it carries the edge scale). Disengaged where fewer than 3 donors reach it.
    std::optional<double> neighbours_median(double frequency) const;
    // log10 f, moved to just inside the gap (1e-6 decades past the edge) when it is not already
    // further: the donors are read as they are beyond the edge, where the join is made.
    double inside_gap(double logFrequency) const;
    const std::string& method(ComplexPermeabilitySource valueSource) const;
};

/**
 * @brief Everything estimated for one material: one entry per part and side its table does not
 * reach. Built once per material (ComplexPermeability caches it) and labelled.
 */
struct ComplexPermeabilityEstimate {
    std::string materialName;
    std::string materialClass;                   // e.g. "FERRITE MN_ZN"
    std::pair<double, double> realSpan;          // tabulated mu', Hz
    std::pair<double, double> imaginarySpan;     // tabulated mu'', Hz
    size_t numberClassDonors = 0;                // class members usable as donors
    std::vector<std::string> snoekScreened;      // class members screened out of the donor pool
    std::vector<ComplexPermeabilityGapEstimate> gaps;

    const ComplexPermeabilityGapEstimate& gap(ComplexPermeabilityPart part, ComplexPermeabilityGapSide side) const;
    std::string summary() const;  // the one-per-material WARNING
};

/**
 * @brief A stretch of frequency in which a part of the complex permeability is estimated.
 * fromFrequency is 0 and toFrequency +infinity for the open ends.
 */
struct ComplexPermeabilityEstimatedRange {
    double fromFrequency;
    double toFrequency;
    ComplexPermeabilityPart part;
    ComplexPermeabilitySource source;
    std::string materialClass;
    std::vector<std::string> donors;
    std::string method;
};

/**
 * @brief Builds the estimate of every part and side that the material's MAS complex-permeability
 * table does not reach, from the other materials of its class in the loaded core-material database.
 * Throws ComplexPermeabilityEstimatorRefusal when it cannot (see above), MaterialException for
 * malformed data (a tabulated frequency that is not a finite positive number).
 */
ComplexPermeabilityEstimate estimate_complex_permeability(const CoreMaterial& coreMaterial);

std::vector<ComplexPermeabilityEstimatedRange> get_estimated_ranges(const ComplexPermeabilityEstimate& estimate);

}  // namespace OpenMagnetics
