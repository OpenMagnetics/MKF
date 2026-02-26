#pragma once
#include "Constants.h"
#include "Defaults.h"
#include "Models.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Wire.h"

#include "processors/Inputs.h"
#include "support/CoilMesher.h"
#include <MAS.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <streambuf>
#include <vector>
#include "support/Exceptions.h"

using namespace MAS; // QUAL-001 TODO: qualify types and remove

namespace OpenMagnetics {

class MagneticFieldStrengthModel {
    public:
        std::vector<Wire> _wirePerWinding;
        virtual ComplexFieldPoint get_magnetic_field_strength_between_two_points(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, std::optional<size_t> inducingWireIndex = std::nullopt) = 0;
};

class MagneticFieldStrengthFringingEffectModel {
    public:
        virtual FieldPoint get_equivalent_inducing_point_for_gap(CoreGap gap, double magneticFieldStrengthGap) = 0;
        virtual ComplexFieldPoint get_magnetic_field_strength_between_gap_and_point(CoreGap gap, double magneticFieldStrengthGap, FieldPoint inducedFieldPoint) = 0;
};

class MagneticField {
    private:
    protected:
        MagneticFieldStrengthModels _magneticFieldStrengthModel;
        MagneticFieldStrengthFringingEffectModels _magneticFieldStrengthFringingEffectModel;
        std::shared_ptr<MagneticFieldStrengthModel>  _model;
        std::shared_ptr<MagneticFieldStrengthFringingEffectModel>  _fringingEffectModel;
    public:
        std::vector<Wire> _wirePerWinding;

        MagneticField(MagneticFieldStrengthModels magneticFieldStrengthModel = Defaults().magneticFieldStrengthModelDefault,
                      MagneticFieldStrengthFringingEffectModels magneticFieldStrengthFringingEffectModel = Defaults().magneticFieldStrengthFringingEffectModelDefault) {
            _magneticFieldStrengthModel = magneticFieldStrengthModel;
            _magneticFieldStrengthFringingEffectModel = magneticFieldStrengthFringingEffectModel;
            _model = factory(_magneticFieldStrengthModel);
            _fringingEffectModel = factory(_magneticFieldStrengthFringingEffectModel);
        }

        static SignalDescriptor calculate_magnetic_flux(SignalDescriptor magnetizingCurrent,
                                                        double MagneticField,
                                                        double numberTurns);
        static SignalDescriptor calculate_magnetic_flux_density(SignalDescriptor magneticFlux,
                                                                double area);
        static SignalDescriptor calculate_magnetic_field_strength(SignalDescriptor magneticFluxDensity,
                                                                    double initialPermeability);

        WindingWindowMagneticStrengthFieldOutput calculate_magnetic_field_strength_field(OperatingPoint operatingPoint, Magnetic magnetic, std::optional<Field> externalInducedField = std::nullopt, std::optional<std::vector<int8_t>> customCurrentDirectionPerWinding = std::nullopt, std::optional<CoilMesherModels> coilMesherModel = std::nullopt);

        static std::shared_ptr<MagneticFieldStrengthFringingEffectModel> factory(MagneticFieldStrengthFringingEffectModels modelName);
        static std::shared_ptr<MagneticFieldStrengthModel> factory(MagneticFieldStrengthModels modelName);
        static std::shared_ptr<MagneticFieldStrengthModel> factory();
};



// // Based on Effects of eddy currents in transformer windings by P. L. Dowell.
// // https://sci-hub.st/10.1049/piee.1966.0236
// class MagneticFieldStrengthDowellModel : public MagneticFieldStrengthModel {
//     public:
//         std::string methodName = "Dowell";
//         ComplexFieldPoint get_magnetic_field_strength_between_two_points(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, std::optional<size_t> inducingWireIndex = std::nullopt);
// };



// Based on Improved Analytical Calculation of High Frequency Winding Losses in Planar Inductors by Xiaohui Wang
// https://sci-hub.st/10.1109/ECCE.2018.8558397
class MagneticFieldStrengthWangModel : public MagneticFieldStrengthModel {
    public:
        std::string methodName = "Wang";
        ComplexFieldPoint get_magnetic_field_strength_between_two_points(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, std::optional<size_t> inducingWireIndex = std::nullopt);
};


// Based on Analysis and Computation of Electric and Magnetic Field Problems By Kenneth John Binns and P. J. Lawrenson
// https://library.lol/main/78A1F466FA9F3EFBCB6165283FC346B6
// Equation 3.34 and 5.4
class MagneticFieldStrengthBinnsLawrensonModel : public MagneticFieldStrengthModel {
    public:
        std::string methodName = "BinnsLawrenson";
        ComplexFieldPoint get_magnetic_field_strength_between_two_points(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, std::optional<size_t> inducingWireIndex = std::nullopt);
};


// Based on Eddy currents by Jiří Lammeraner
// https://archive.org/details/eddycurrents0000lamm
class MagneticFieldStrengthLammeranerModel : public MagneticFieldStrengthModel {
    public:
        std::string methodName = "Lammeraner";
        ComplexFieldPoint get_magnetic_field_strength_between_two_points(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, std::optional<size_t> inducingWireIndex = std::nullopt);
};


// Based on Fringing Field Formulas and Winding Loss Due to an Air Gap by Waseem A. Roshen
// https://sci-hub.st/10.1109/tmag.2007.898908
// https://sci-hub.st/10.1109/tmag.2008.2002302
class MagneticFieldStrengthRoshenModel : public MagneticFieldStrengthFringingEffectModel {
    public:
        std::string methodName = "Roshen";
        ComplexFieldPoint get_magnetic_field_strength_between_gap_and_point(CoreGap gap, double magneticFieldStrengthGap, FieldPoint inducedFieldPoint);
        FieldPoint get_equivalent_inducing_point_for_gap([[maybe_unused]]CoreGap gap, [[maybe_unused]]double magneticFieldStrengthGap) {
            throw std::runtime_error("Fringing field not implemented for this model");
        }
};


// Based on Induktivitäten in der Leistungselektronik: Spulen, Trafos und ihre parasitären Eigenschaften by Manfred Albach
// https://libgen.rocks/get.php?md5=94b7f2906f53602f19892d7f1dabd929&key=YMKCEJOWB653PYLL
// https://sci-hub.st/10.1109/tpel.2011.2143729
class MagneticFieldStrengthAlbachModel : public MagneticFieldStrengthFringingEffectModel {
    public:
        std::string methodName = "Albach";
        ComplexFieldPoint get_magnetic_field_strength_between_gap_and_point([[maybe_unused]]CoreGap gap, [[maybe_unused]]double magneticFieldStrengthGap, [[maybe_unused]]FieldPoint inducedFieldPoint) {
            throw std::runtime_error("Fringing field not implemented for this model");
        }
        FieldPoint get_equivalent_inducing_point_for_gap(CoreGap gap, double magneticFieldStrengthGap);
};

// ============================================================================
// SULLIVAN Fringing Field Model (2D Image Method / Biot-Savart)
// ============================================================================
/**
 * @brief Fringing field model based on the method of images (Biot-Savart approach)
 *
 * Based on the "shapeopt" MATLAB tool by C.R. Sullivan et al. at Dartmouth College.
 * Reference: http://thayer.dartmouth.edu/inductor
 *
 * Key publications:
 * - J. Hu, C.R. Sullivan, "Optimization of shapes for round-wire high-frequency
 *   gapped-inductor windings", IEEE IAS Annual Meeting, 1998
 * - C.R. Sullivan, "Optimal Choice for Number of Strands in a Litz-Wire
 *   Transformer Winding", IEEE Trans. Power Electron., 14(2), March 1999
 *
 * This model uses the 2D method of images to compute fringing magnetic field
 * strength (H) from an air gap. The gap MMF is modeled as distributed current
 * filaments across the gap width. Each filament is placed at the gap position
 * along with its image (reflected about the centerpost face), and these pairs
 * are replicated across image units of the winding window to enforce boundary
 * conditions imposed by the high-permeability core walls.
 *
 * The field is computed via 2D Biot-Savart superposition:
 *   B = (mu_0 * I) / (2*pi) * sum_images( (R - R_source) / |R - R_source|^2 )
 * converted to H = B / mu_0
 *
 * Unlike Roshen (closed-form) or Albach (equivalent current loop), this model
 * captures the effect of multiple reflections from all core walls, making it
 * more accurate for narrow winding windows or points close to core walls.
 */
class MagneticFieldStrengthSullivanModel : public MagneticFieldStrengthFringingEffectModel {
 public:
    std::string methodName = "Sullivan";

    /**
     * @brief Compute fringing H-field at an arbitrary point due to one gap
     *
     * Uses 2D Biot-Savart with method of images. The gap is subdivided into
     * _gapDivisions filaments, each replicated across (2*_imageUnitsX+1) x
     * (2*_imageUnitsY+1) image cells.
     *
     * @param gap                      CoreGap with coordinates, length, section_dimensions
     * @param magneticFieldStrengthGap Peak H in the gap (A/m), i.e. B_gap / mu_0
     * @param inducedFieldPoint        Point where the field is evaluated
     * @return ComplexFieldPoint with real = Hx (radial), imaginary = Hy (axial)
     */
    ComplexFieldPoint get_magnetic_field_strength_between_gap_and_point(
        CoreGap gap, double magneticFieldStrengthGap, FieldPoint inducedFieldPoint);

    /**
     * @brief Not used - Sullivan computes field directly, not via equivalent point
     */
    FieldPoint get_equivalent_inducing_point_for_gap(
        [[maybe_unused]] CoreGap gap, [[maybe_unused]] double magneticFieldStrengthGap) {
        throw std::runtime_error("Sullivan fringing model computes field directly, not via equivalent point");
    }

 private:
    /// Number of subdivisions of the gap length (higher = more accurate, slower)
    int _gapDivisions = 10;
    /// Number of image unit repetitions in the X direction (perpendicular to gap)
    int _imageUnitsX = 2;
    /// Number of image unit repetitions in the Y direction (parallel to gap)
    int _imageUnitsY = 2;
};



// ============================================================================
// ALBACH H-Field Model (Air Coil / Biot-Savart approach)
// ============================================================================

/**
 * @brief Structure representing a turn's position for the ALBACH H-field model
 * 
 * For round wires: width and height should be 0 (point filament approximation)
 * For rectangular wires: width = radial extent, height = axial extent
 *   The current is assumed uniformly distributed across the cross-section
 *   and is modeled using filamentary subdivision (multiple circular filaments)
 */
struct AlbachTurnPosition {
    double r;           ///< Radial position of conductor center (m)
    double z;           ///< Axial position of conductor center (m)
    double current;     ///< Current amplitude (A)
    size_t turnIndex;   ///< Index in the coil's turn list
    
    // For rectangular wire support (set to 0 for round wires)
    double width = 0;   ///< Radial extent of conductor cross-section (m)
    double height = 0;  ///< Axial extent of conductor cross-section (m)
    
    // For frequency-dependent current distribution
    double skinDepth = 1e9;  ///< Skin depth at current frequency (m), large default = uniform distribution

    /// @brief Check if this turn represents a rectangular conductor
    bool isRectangular() const { return width > 1e-10 && height > 1e-10; }
};

/**
 * @brief Magnetic field strength model using air coil (Biot-Savart) approach
 * 
 * Based on M. Albach, "Two-dimensional calculation of winding losses in transformers",
 * PESC 2000.
 * 
 * This model computes the H-field contribution from TURNS only using the
 * analytical Biot-Savart formula with elliptic integrals for circular current
 * filaments. Gap fringing effects are handled separately by the fringing
 * effect models (ROSHEN or ALBACH fringing).
 */
class MagneticFieldStrengthAlbach2DModel : public MagneticFieldStrengthModel {
public:
    std::string methodName = "Albach";
    
    /**
     * @brief Calculate H field between inducing and induced points
     * 
     * This model calculates field from all turns at once via calculateTotalFieldAtPoint().
     * This per-turn-pair method should not be called and will throw an error.
     */
    ComplexFieldPoint get_magnetic_field_strength_between_two_points(
        FieldPoint inducingFieldPoint, 
        FieldPoint inducedFieldPoint, 
        std::optional<size_t> inducingWireIndex = std::nullopt
    ) override;
    
    /**
     * @brief Set the turn positions for field calculation
     */
    void setTurns(const std::vector<AlbachTurnPosition>& turns) {
        _turns = turns;
    }
    
    /**
     * @brief Update turn currents for a specific harmonic
     * This allows reusing the setup with different current values
     */
    void updateTurnCurrents(const std::vector<double>& currents) {
        for (size_t i = 0; i < _turns.size() && i < currents.size(); ++i) {
            _turns[i].current = currents[i];
        }
    }

    /**
     * @brief Update skin depths for all turns based on frequency
     * 
     * At high frequency, current concentrates at conductor edges due to skin effect.
     * 
     * @param skinDepth Skin depth at current frequency (m)
     */
    void updateSkinDepths(double skinDepth) {
        for (auto& turn : _turns) {
            turn.skinDepth = skinDepth;
        }
    }
    
    /**
     * @brief Calculate total field at an induced point from all turns
     * This is more efficient than the per-turn-pair interface
     */
    ComplexFieldPoint calculateTotalFieldAtPoint(FieldPoint inducedFieldPoint);
    
    /**
     * @brief Setup the model from a Magnetic component
     * 
     * Extracts turn positions from the coil and sets up for field calculation.
     * 
     * @param magnetic The magnetic component (core + coil)
     * @param wirePerWinding Wire definitions for each winding
     */
    void setupFromMagnetic(Magnetic magnetic, const std::vector<Wire>& wirePerWinding);
    
private:
    std::vector<AlbachTurnPosition> _turns;
    
    /**
     * @brief Calculate H field at point (r, z) from all turns using Biot-Savart
     * 
     * @param r Radial coordinate of field point (m)
     * @param z Axial coordinate of field point (m)
     * @return pair<H_r, H_z> Magnetic field components in A/m
     */
    std::pair<double, double> calculateMagneticField(double r, double z) const;
};

} // namespace OpenMagnetics