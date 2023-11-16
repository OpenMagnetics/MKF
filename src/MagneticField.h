#pragma once
#include "Constants.h"
#include "Defaults.h"
#include "Models.h"
#include "MagneticWrapper.h"

#include <InputsWrapper.h>
#include <MAS.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <map>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {

class MagneticFieldStrengthModel {
  public:
    virtual ComplexFieldPoint get_magnetic_field_strength_between_two_points(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, std::optional<WireWrapper> inducingWire = std::nullopt) = 0;
    virtual ComplexFieldPoint get_magnetic_field_strength_between_gap_and_point(CoreGap gap, int64_t numberTurns, double magnetizingCurrentPeak, FieldPoint inducedFieldPoint) = 0;
};

class MagneticField {
    private:
    protected:
        MagneticFieldStrengthModels _magneticFieldStrengthModel;
        std::shared_ptr<MagneticFieldStrengthModel>  _model;
        std::shared_ptr<MagneticFieldStrengthModel>  _fringingEffectModel;
        bool _includeFringing = true;
        int _mirroringDimension = 0;
    public:

        MagneticField(MagneticFieldStrengthModels magneticFieldStrengthModel = Defaults().magneticFieldStrengthModelDefault) {
            _magneticFieldStrengthModel = magneticFieldStrengthModel;
            _model = factory(_magneticFieldStrengthModel);
            _fringingEffectModel = factory(MagneticFieldStrengthModels::ROSHEN);
        }

        void set_fringing_effect(bool value) {
            _includeFringing = value;
        }

        void set_mirroring_dimension(int mirroringDimension) {
            _mirroringDimension = mirroringDimension;
        }

        static SignalDescriptor calculate_magnetic_flux(SignalDescriptor magnetizingCurrent,
                                                        double MagneticField,
                                                        double numberTurns,
                                                        double frequency);
        static SignalDescriptor calculate_magnetic_flux_density(SignalDescriptor magneticFlux,
                                                                double area,
                                                                double frequency);
        static SignalDescriptor calculate_magnetic_field_strength(SignalDescriptor magneticFluxDensity,
                                                                    double initialPermeability,
                                                                    double frequency);

        WindingWindowMagneticStrengthFieldOutput calculate_magnetic_field_strength_field(OperatingPoint operatingPoint, MagneticWrapper magnetic, std::optional<std::vector<FieldPoint>> externalInducedFieldPoints = std::nullopt);

        static std::shared_ptr<MagneticFieldStrengthModel> factory(MagneticFieldStrengthModels modelName);
        static std::shared_ptr<MagneticFieldStrengthModel> factory();
};



// Based on Effects of eddy currents in transformer windings by P. L. Dowell.
// https://sci-hub.wf/10.1049/piee.1966.0236
class MagneticFieldStrengthDowellModel : public MagneticFieldStrengthModel {
  public:
    ComplexFieldPoint get_magnetic_field_strength_between_two_points(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, std::optional<WireWrapper> inducingWire = std::nullopt);
    ComplexFieldPoint get_magnetic_field_strength_between_gap_and_point([[maybe_unused]]CoreGap gap, [[maybe_unused]]int64_t numberTurns, [[maybe_unused]]double magnetizingCurrentPeak, [[maybe_unused]] FieldPoint inducedFieldPoint) {
        throw std::runtime_error("Fringing field not implemented for this model");
    }
};


// Based on Analysis and Computation of Electric and Magnetic Field Problems By Kenneth John Binns and P. J. Lawrenson
// https://library.lol/main/78A1F466FA9F3EFBCB6165283FC346B6
class MagneticFieldStrengthBinnsLawrensonModel : public MagneticFieldStrengthModel {
  public:
    ComplexFieldPoint get_magnetic_field_strength_between_two_points(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, std::optional<WireWrapper> inducingWire = std::nullopt);
    ComplexFieldPoint get_magnetic_field_strength_between_gap_and_point([[maybe_unused]]CoreGap gap, [[maybe_unused]]int64_t numberTurns, [[maybe_unused]]double magnetizingCurrentPeak, [[maybe_unused]] FieldPoint inducedFieldPoint) {
        throw std::runtime_error("Fringing field not implemented for this model");
    }
};


// Based on Eddy currents by Jiří Lammeraner
// https://archive.org/details/eddycurrents0000lamm
class MagneticFieldStrengthLammeranerModel : public MagneticFieldStrengthModel {
  public:
    ComplexFieldPoint get_magnetic_field_strength_between_two_points(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, std::optional<WireWrapper> inducingWire = std::nullopt);
    ComplexFieldPoint get_magnetic_field_strength_between_gap_and_point([[maybe_unused]]CoreGap gap, [[maybe_unused]]int64_t numberTurns, [[maybe_unused]]double magnetizingCurrentPeak, [[maybe_unused]] FieldPoint inducedFieldPoint) {
        throw std::runtime_error("Fringing field not implemented for this model");
    }
};


// Based on Fringing Field Formulas and Winding Loss Due to an Air Gap by Waseem A. Roshen
// https://sci-hub.wf/10.1109/tmag.2007.898908
// https://sci-hub.st/10.1109/tmag.2008.2002302
class MagneticFieldStrengthRoshenModel : public MagneticFieldStrengthModel {
  public:
    ComplexFieldPoint get_magnetic_field_strength_between_two_points([[maybe_unused]]FieldPoint inducingFieldPoint, [[maybe_unused]]FieldPoint inducedFieldPoint, [[maybe_unused]]std::optional<WireWrapper> inducingWire = std::nullopt) {
        throw std::runtime_error("Point to point not implemented for this model");
    }
    ComplexFieldPoint get_magnetic_field_strength_between_gap_and_point(CoreGap gap, int64_t numberTurns, double magnetizingCurrentPeak, FieldPoint inducedFieldPoint);
};

} // namespace OpenMagnetics