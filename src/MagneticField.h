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
        bool _includeFringing = true;
        int _mirroringDimension = 0;
        double _windingLossesHarmonicAmplitudeThreshold = Defaults().windingLossesHarmonicAmplitudeThreshold;
    public:

        MagneticField(MagneticFieldStrengthModels magneticFieldStrengthModel = Defaults().magneticFieldStrengthModelDefault,
                      MagneticFieldStrengthFringingEffectModels magneticFieldStrengthFringingEffectModel = Defaults().magneticFieldStrengthFringingEffectModelDefault) {
            _magneticFieldStrengthModel = magneticFieldStrengthModel;
            _magneticFieldStrengthFringingEffectModel = magneticFieldStrengthFringingEffectModel;
            _model = factory(_magneticFieldStrengthModel);
            _fringingEffectModel = factory(_magneticFieldStrengthFringingEffectModel);
        }

        void set_fringing_effect(bool value) {
            _includeFringing = value;
        }

        void set_mirroring_dimension(int mirroringDimension) {
            _mirroringDimension = mirroringDimension;
        }

        void set_winding_losses_harmonic_amplitude_threshold(double windingLossesHarmonicAmplitudeThreshold) {
            _windingLossesHarmonicAmplitudeThreshold = windingLossesHarmonicAmplitudeThreshold;
        }

        static SignalDescriptor calculate_magnetic_flux(SignalDescriptor magnetizingCurrent,
                                                        double MagneticField,
                                                        double numberTurns);
        static SignalDescriptor calculate_magnetic_flux_density(SignalDescriptor magneticFlux,
                                                                double area);
        static SignalDescriptor calculate_magnetic_field_strength(SignalDescriptor magneticFluxDensity,
                                                                    double initialPermeability);

        WindingWindowMagneticStrengthFieldOutput calculate_magnetic_field_strength_field(OperatingPoint operatingPoint, MagneticWrapper magnetic, std::optional<Field> externalInducedField = std::nullopt);

        static std::shared_ptr<MagneticFieldStrengthFringingEffectModel> factory(MagneticFieldStrengthFringingEffectModels modelName);
        static std::shared_ptr<MagneticFieldStrengthModel> factory(MagneticFieldStrengthModels modelName);
        static std::shared_ptr<MagneticFieldStrengthModel> factory();
};



// Based on Effects of eddy currents in transformer windings by P. L. Dowell.
// https://sci-hub.wf/10.1049/piee.1966.0236
class MagneticFieldStrengthDowellModel : public MagneticFieldStrengthModel {
    public:
        ComplexFieldPoint get_magnetic_field_strength_between_two_points(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, std::optional<WireWrapper> inducingWire = std::nullopt);
};


// Based on Analysis and Computation of Electric and Magnetic Field Problems By Kenneth John Binns and P. J. Lawrenson
// https://library.lol/main/78A1F466FA9F3EFBCB6165283FC346B6
class MagneticFieldStrengthBinnsLawrensonModel : public MagneticFieldStrengthModel {
    public:
        ComplexFieldPoint get_magnetic_field_strength_between_two_points(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, std::optional<WireWrapper> inducingWire = std::nullopt);
};


// Based on Eddy currents by Jiří Lammeraner
// https://archive.org/details/eddycurrents0000lamm
class MagneticFieldStrengthLammeranerModel : public MagneticFieldStrengthModel {
    public:
        ComplexFieldPoint get_magnetic_field_strength_between_two_points(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, std::optional<WireWrapper> inducingWire = std::nullopt);
};


// Based on Fringing Field Formulas and Winding Loss Due to an Air Gap by Waseem A. Roshen
// https://sci-hub.wf/10.1109/tmag.2007.898908
// https://sci-hub.st/10.1109/tmag.2008.2002302
class MagneticFieldStrengthRoshenModel : public MagneticFieldStrengthFringingEffectModel {
    public:
        ComplexFieldPoint get_magnetic_field_strength_between_gap_and_point(CoreGap gap, double magneticFieldStrengthGap, FieldPoint inducedFieldPoint);
        FieldPoint get_equivalent_inducing_point_for_gap([[maybe_unused]]CoreGap gap, [[maybe_unused]]double magneticFieldStrengthGap) {
            throw std::runtime_error("Fringing field not implemented for this model");
        }
};


// Based on Induktivitäten in der Leistungselektronik: Spulen, Trafos und ihre parasitären Eigenschaften by Manfred Albach
// https://libgen.rocks/get.php?md5=94b7f2906f53602f19892d7f1dabd929&key=YMKCEJOWB653PYLL
// https://sci-hub.wf/10.1109/tpel.2011.2143729
class MagneticFieldStrengthAlbachModel : public MagneticFieldStrengthFringingEffectModel {
    public:
        ComplexFieldPoint get_magnetic_field_strength_between_gap_and_point([[maybe_unused]]CoreGap gap, [[maybe_unused]]double magneticFieldStrengthGap, [[maybe_unused]]FieldPoint inducedFieldPoint) {
            throw std::runtime_error("Fringing field not implemented for this model");
        }
        FieldPoint get_equivalent_inducing_point_for_gap(CoreGap gap, double magneticFieldStrengthGap);
};


} // namespace OpenMagnetics