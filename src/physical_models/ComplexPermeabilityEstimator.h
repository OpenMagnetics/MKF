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
 *   its mu' and mu'' are scaled, to fit the target's OWN data on the known side of the gap, weighted
 *   down with distance from the edge (exp(-d / 1 decade)), plus the other part inside the gap, where
 *   that part is tabulated. The estimate at a frequency is the median of the 5 best-matching donors
 *   whose shifted data reach it (at least 3 required), scaled so that it meets the material's own
 *   last tabulated value at the edge exactly (the donors give the shape, the material the level).
 * - Beyond the donors' reach, the estimate continues from its outermost point along the class's tail
 *   slope, learned from the catalogue (median log-log slope of the donors' outermost decade).
 * - Below its data, mu' is held at its lowest tabulated value (the catalogue shows mu' flat there).
 * - Donors whose mu_i * f(mu'' peak) is more than 2 decades off their class median are screened out
 *   (19 DMEGC records and ACME N5 sit ~3 decades low: frequencies that look like kHz stored as Hz).
 * Hold-out errors of |Z| (toroid) with data removed above 1 MHz, MnZn: median 17 %, P90 73 % (C++ port,
 * with the edge scale and the distance weights; the proposal's Python: 19 % / 68 %), against
 * 182 % / 1696 % for the best simple baseline.
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
    HELD_FLAT              // mu' below its data: its lowest tabulated value
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
    double reachValue = 0;                       // the estimate there
    double tailSlope = 0;                        // learned: d log|value| / d log f past the reach
    size_t tailSlopeDonors = 0;                  // donors the tail slope is the median of
    // D7 (Alf, 2026-10-14): every donor's gain carries edgeScale = the part's own last tabulated value /
    // the donors' median there, so the estimate meets the material's own data at the edge exactly.
    double edgeScale = 1;
    // edgeScale < 0: the donors' median at the edge has the opposite sign of the own mu' there, and the
    // multiplicative scale flips the estimated mu' (reported to Alf, ABT #1861; not hidden).
    bool edgeSignMismatch = false;
    std::string neighboursMethod;                // human-readable, as CoreHeatConductivityEstimate::method
    std::string tailMethod;
    // Non-empty when this part and side cannot be estimated (e.g. no donor's data reach past the
    // edge): evaluate() throws ComplexPermeabilityEstimatorRefusal with it. The other gaps stand.
    std::string refusal;

    // The estimated value at a frequency beyond the edge, and where it came from.
    std::pair<double, ComplexPermeabilitySource> evaluate(double frequency) const;
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
