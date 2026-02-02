#include "support/Utils.h"
#include "support/Settings.h"

namespace OpenMagnetics {

    Settings& Settings::GetInstance()
    {
        static Settings instance;
        return instance;
    }

    Settings::Settings() {
        _inputsNumberPointsSampledWaveforms = constants.numberPointsSampledWaveforms;
        _painterMirroringDimension = defaults.magneticFieldMirroringDimension;
        _magneticFieldMirroringDimension = defaults.magneticFieldMirroringDimension;
        _harmonicAmplitudeThreshold = defaults.harmonicAmplitudeThreshold;
        _coreLossesModelNames = {defaults.coreLossesModelDefault, CoreLossesModels::PROPRIETARY, CoreLossesModels::LOSS_FACTOR, CoreLossesModels::STEINMETZ, CoreLossesModels::ROSHEN};
        _coreAdviserMaximumMagneticsAfterFiltering = defaults.coreAdviserMaximumMagneticsAfterFiltering;
    }


    void Settings::reset() {

        _useToroidalCores = true;
        _useConcentricCores = true;

        _inputsTrimHarmonics = true;

        _inputsNumberPointsSampledWaveforms = Constants().numberPointsSampledWaveforms;

        _magnetizingInductanceIncludeAirInductance = false;

        _coilAllowMarginTape = true;
        _coilAllowInsulatedWire = true;
        _coilFillSectionsWithMarginTape = false;
        _coilWindEvenIfNotFit = false;
        _coilDelimitAndCompact = true;
        _coilTryRewind = true;
        _coilIncludeAdditionalCoordinates = true;
        _coilEqualizeMargins = true;
        _coilOnlyOneTurnPerLayerInContiguousRectangular = false;
        _coilMaximumLayersPlanar = 32;

        _useOnlyCoresInStock = true;

        _painterNumberPointsX = 25;
        _painterNumberPointsY = 50;
        _painterMirroringDimension = Defaults().magneticFieldMirroringDimension;
        _painterMode = PainterModes::CONTOUR;
        _painterLogarithmicScale = false;
        _painterIncludeFringing = true;
        _painterDrawSpacer = true;
        _painterSimpleLitz = true;
        _painterAdvancedLitz = false;
        _painterMaximumValueColorbar = std::nullopt;
        _painterMinimumValueColorbar = std::nullopt;
        _painterColorFerrite = "0x7b7c7d";
        _painterColorBobbin = "0x539796";
        _painterColorCopper = "0xb87333";
        _painterColorInsulation = "0xfff05b";
        _painterColorFr4 = "0x008000";
        _painterColorMargin = "0xfff05b";
        _painterColorEnamel = "0xc63032";
        _painterColorFEP = "0x252525";
        _painterColorETFE = "0xb42811";
        _painterColorTCA = "0x696969";
        _painterColorTCA = "0x696969";
        _painterColorSilk = "0xe7e7e8";
        _painterColorSpacer = "0x3b3b3b";
        _painterColorLines = "0x010000";
        _painterColorText = "0x000000";
        _painterColorCurrentDensity = "0x0892D0";
        _painterCciCoordinatesPath = std::string{selfFilePath}.substr(0, std::string{selfFilePath}.rfind("/")).append("/../../cci_coords/coordinates/");
        _painterColorMagneticFieldMinimum = "0x2b35f5";
        _painterColorMagneticFieldMaximum = "0xe84922";

        _magneticFieldNumberPointsX = 25;
        _magneticFieldNumberPointsY = 50;
        _magneticFieldMirroringDimension = Defaults().magneticFieldMirroringDimension;
        _magneticFieldIncludeFringing = true;

        _coilMesherInsideTurnsFactor = 1.05;

        _leakageInductanceGridAutoScaling = true;
        _leakageInductanceGridPrecisionLevelPlanar = 3;
        _leakageInductanceGridPrecisionLevelWound = 1;
        
        _coilAdviserMaximumNumberWires = 100;
        _coreAdviserIncludeStacks = true;
        _coreAdviserIncludeDistributedGaps = true;
        _coreAdviserIncludeMargin = false;
        _coreAdviserEnableIntermediatePruning = true;
        _coreAdviserMaximumMagneticsAfterFiltering = Defaults().coreAdviserMaximumMagneticsAfterFiltering;


        _wireAdviserIncludePlanar = false;
        _wireAdviserIncludeFoil = false;
        _wireAdviserIncludeRectangular = true;
        _wireAdviserIncludeLitz = true;
        _wireAdviserIncludeRound = true;
        _wireAdviserAllowRectangularInToroidalCores = false;

        _harmonicAmplitudeThresholdQuickMode = true;
        _harmonicAmplitudeThreshold = Defaults().harmonicAmplitudeThreshold;

        _verbose = false;

        _preferredCoreMaterialFerriteManufacturer = "Fair-Rite";
        _preferredCoreMaterialPowderManufacturer = "Micrometals";

        _coreCrossReferencerAllowDifferentCoreMaterialType = false;
    }

    bool Settings::get_verbose() const {
        return _verbose;
    }
    void Settings::set_verbose(bool value) {
        _verbose = value;
    }

    bool Settings::get_use_toroidal_cores() const {
        return _useToroidalCores;
    }
    void Settings::set_use_toroidal_cores(bool value) {
        _useToroidalCores = value;
    }

    bool Settings::get_use_concentric_cores() const {
        return _useConcentricCores;
    }
    void Settings::set_use_concentric_cores(bool value) {
        _useConcentricCores = value;
    }

    bool Settings::get_inputs_trim_harmonics() const {
        return _inputsTrimHarmonics;
    }
    void Settings::set_inputs_trim_harmonics(bool value) {
        _inputsTrimHarmonics = value;
    }

    size_t Settings::get_inputs_number_points_sampled_waveforms() const {
        return _inputsNumberPointsSampledWaveforms;
    }
    void Settings::set_inputs_number_points_sampled_waveforms(size_t value) {
        _inputsNumberPointsSampledWaveforms = value;
    }

    bool Settings::get_magnetizing_inductance_include_air_inductance() const {
        return _magnetizingInductanceIncludeAirInductance;
    }
    void Settings::set_magnetizing_inductance_include_air_inductance(bool value) {
        _magnetizingInductanceIncludeAirInductance = value;
    }

    bool Settings::get_coil_allow_margin_tape() const {
        return _coilAllowMarginTape;
    }
    void Settings::set_coil_allow_margin_tape(bool value) {
        _coilAllowMarginTape = value;
    }

    bool Settings::get_coil_allow_insulated_wire() const {
        return _coilAllowInsulatedWire;
    }
    void Settings::set_coil_allow_insulated_wire(bool value) {
        _coilAllowInsulatedWire = value;
    }

    bool Settings::get_coil_fill_sections_with_margin_tape() const {
        return _coilFillSectionsWithMarginTape;
    }
    void Settings::set_coil_fill_sections_with_margin_tape(bool value) {
        _coilFillSectionsWithMarginTape = value;
    }

    bool Settings::get_coil_wind_even_if_not_fit() const {
        return _coilWindEvenIfNotFit;
    }
    void Settings::set_coil_wind_even_if_not_fit(bool value) {
        _coilWindEvenIfNotFit = value;
    }

    bool Settings::get_coil_delimit_and_compact() const {
        return _coilDelimitAndCompact;
    }
    void Settings::set_coil_delimit_and_compact(bool value) {
        _coilDelimitAndCompact = value;
    }

    bool Settings::get_coil_try_rewind() const {
        return _coilTryRewind;
    }
    void Settings::set_coil_try_rewind(bool value) {
        _coilTryRewind = value;
    }

    bool Settings::get_coil_include_additional_coordinates() const {
        return _coilIncludeAdditionalCoordinates;
    }
    void Settings::set_coil_include_additional_coordinates(bool value) {
        _coilIncludeAdditionalCoordinates = value;
    }

    bool Settings::get_coil_equalize_margins() const {
        return _coilEqualizeMargins;
    }
    void Settings::set_coil_equalize_margins(bool value) {
        _coilEqualizeMargins = value;
    }

    bool Settings::get_coil_only_one_turn_per_layer_in_contiguous_rectangular() const {
        return _coilOnlyOneTurnPerLayerInContiguousRectangular;
    }
    void Settings::set_coil_only_one_turn_per_layer_in_contiguous_rectangular(bool value) {
        _coilOnlyOneTurnPerLayerInContiguousRectangular = value;
    }

    size_t Settings::get_coil_maximum_layers_planar() const {
        return _coilMaximumLayersPlanar;
    }
    void Settings::set_coil_maximum_layers_planar(size_t value) {
        _coilMaximumLayersPlanar = value;
    }

    bool Settings::get_use_only_cores_in_stock() const {
        return _useOnlyCoresInStock;
    }
    void Settings::set_use_only_cores_in_stock(bool value) {
        _useOnlyCoresInStock = value;
    }

    size_t Settings::get_painter_number_points_x() const {
        return _painterNumberPointsX;
    }
    void Settings::set_painter_number_points_x(size_t value) {
        _painterNumberPointsX = value;
    }

    size_t Settings::get_painter_number_points_y() const {
        return _painterNumberPointsY;
    }
    void Settings::set_painter_number_points_y(size_t value) {
        _painterNumberPointsY = value;
    }

    PainterModes Settings::get_painter_mode() const {
        return _painterMode;
    }
    void Settings::set_painter_mode(PainterModes value) {
        _painterMode = value;
    }

    bool Settings::get_painter_logarithmic_scale() const {
        return _painterLogarithmicScale;
    }
    void Settings::set_painter_logarithmic_scale(bool value) {
        _painterLogarithmicScale = value;
    }

    bool Settings::get_painter_include_fringing() const {
        return _painterIncludeFringing;
    }
    void Settings::set_painter_include_fringing(bool value) {
        _painterIncludeFringing = value;
    }

    bool Settings::get_painter_draw_spacer() const {
        return _painterDrawSpacer;
    }
    void Settings::set_painter_draw_spacer(bool value) {
        _painterDrawSpacer = value;
    }

    bool Settings::get_painter_simple_litz() const {
        return _painterSimpleLitz;
    }
    void Settings::set_painter_simple_litz(bool value) {
        _painterSimpleLitz = value;
    }

    bool Settings::get_painter_advanced_litz() const {
        return _painterAdvancedLitz;
    }
    void Settings::set_painter_advanced_litz(bool value) {
        _painterAdvancedLitz = value;
    }

    std::optional<double> Settings::get_painter_maximum_value_colorbar() const {
        return _painterMaximumValueColorbar;
    }
    void Settings::set_painter_maximum_value_colorbar(std::optional<double> value) {
        _painterMaximumValueColorbar = value;
    }

    std::optional<double> Settings::get_painter_minimum_value_colorbar() const {
        return _painterMinimumValueColorbar;
    }
    void Settings::set_painter_minimum_value_colorbar(std::optional<double> value) {
        _painterMinimumValueColorbar = value;
    }

    std::string Settings::get_painter_color_ferrite() const {
        return _painterColorFerrite;
    }
    void Settings::set_painter_color_ferrite(std::string value) {
        _painterColorFerrite = value;
    }

    std::string Settings::get_painter_color_bobbin() const {
        return _painterColorBobbin;
    }
    void Settings::set_painter_color_bobbin(std::string value) {
        _painterColorBobbin = value;
    }

    std::string Settings::get_painter_color_copper() const {
        return _painterColorCopper;
    }
    void Settings::set_painter_color_copper(std::string value) {
        _painterColorCopper = value;
    }

    std::string Settings::get_painter_color_insulation() const {
        return _painterColorInsulation;
    }
    void Settings::set_painter_color_insulation(std::string value) {
        _painterColorInsulation = value;
    }

    std::string Settings::get_painter_color_fr4() const {
        return _painterColorFr4;
    }
    void Settings::set_painter_color_fr4(std::string value) {
        _painterColorFr4 = value;
    }

    std::string Settings::get_painter_color_margin() const {
        return _painterColorMargin;
    }
    void Settings::set_painter_color_margin(std::string value) {
        _painterColorMargin = value;
    }

    std::string Settings::get_painter_color_spacer() const {
        return _painterColorSpacer;
    }
    void Settings::set_painter_color_spacer(std::string value) {
        _painterColorSpacer = value;
    }

    std::string Settings::get_painter_color_lines() const {
        return _painterColorLines;
    }
    void Settings::set_painter_color_lines(std::string value) {
        _painterColorLines = value;
    }

    std::string Settings::get_painter_color_text() const {
        return _painterColorText;
    }
    void Settings::set_painter_color_text(std::string value) {
        _painterColorText = value;
    }

    std::string Settings::get_painter_cci_coordinates_path() const {
        return _painterCciCoordinatesPath;
    }
    void Settings::set_painter_cci_coordinates_path(std::string value) {
        _painterCciCoordinatesPath = value;
    }

    std::string Settings::get_painter_color_enamel() const {
        return _painterColorEnamel;
    }

    void Settings::set_painter_color_enamel(std::string value) {
        _painterColorEnamel = value;
    }

    std::string Settings::get_painter_color_fep() const {
        return _painterColorFEP;
    }

    void Settings::set_painter_color_fep(std::string value) {
        _painterColorFEP = value;
    }

    std::string Settings::get_painter_color_etfe() const {
        return _painterColorETFE;
    }

    void Settings::set_painter_color_etfe(std::string value) {
        _painterColorETFE = value;
    }

    std::string Settings::get_painter_color_tca() const {
        return _painterColorTCA;
    }

    void Settings::set_painter_color_tca(std::string value) {
        _painterColorTCA = value;
    }

    std::string Settings::get_painter_color_pfa() const {
        return _painterColorPFA;
    }

    void Settings::set_painter_color_pfa(std::string value) {
        _painterColorPFA = value;
    }

    std::string Settings::get_painter_color_silk() const {
        return _painterColorSilk;
    }

    void Settings::set_painter_color_silk(std::string value) {
        _painterColorSilk = value;
    }

    std::string Settings::get_painter_color_current_density() const {
        return _painterColorCurrentDensity;
    }

    void Settings::set_painter_color_current_density(std::string value) {
        _painterColorCurrentDensity = value;
    }

    std::string Settings::get_painter_color_magnetic_field_minimum() const {
        return _painterColorMagneticFieldMinimum;
    }

    void Settings::set_painter_color_magnetic_field_minimum(std::string value) {
        _painterColorMagneticFieldMinimum = value;
    }

    std::string Settings::get_painter_color_magnetic_field_maximum() const {
        return _painterColorMagneticFieldMaximum;
    }

    void Settings::set_painter_color_magnetic_field_maximum(std::string value) {
        _painterColorMagneticFieldMaximum = value;
    }

    int Settings::get_painter_mirroring_dimension() const {
        return _painterMirroringDimension;
    }
    void Settings::set_painter_mirroring_dimension(int value) {
        _painterMirroringDimension = value;
    }

    size_t Settings::get_magnetic_field_number_points_x() const {
        return _magneticFieldNumberPointsX;
    }
    void Settings::set_magnetic_field_number_points_x(size_t value) {
        _magneticFieldNumberPointsX = value;
    }

    size_t Settings::get_magnetic_field_number_points_y() const {
        return _magneticFieldNumberPointsY;
    }
    void Settings::set_magnetic_field_number_points_y(size_t value) {
        _magneticFieldNumberPointsY = value;
    }

    int Settings::get_magnetic_field_mirroring_dimension() const {
        return _magneticFieldMirroringDimension;
    }
    void Settings::set_magnetic_field_mirroring_dimension(int value) {
        _magneticFieldMirroringDimension = value;
    }

    bool Settings::get_magnetic_field_include_fringing() const {
        return _magneticFieldIncludeFringing;
    }
    void Settings::set_magnetic_field_include_fringing(bool value) {
        _magneticFieldIncludeFringing = value;
    }

    double Settings::get_coil_mesher_inside_turns_factor() const {
        return _coilMesherInsideTurnsFactor;
    }
    void Settings::set_coil_mesher_inside_turns_factor(double value) {
        _coilMesherInsideTurnsFactor = value;
    }

    bool Settings::get_leakage_inductance_grid_auto_scaling() const {
        return _leakageInductanceGridAutoScaling;
    }
    void Settings::set_leakage_inductance_grid_auto_scaling(bool value) {
        _leakageInductanceGridAutoScaling = value;
    }

    double Settings::get_leakage_inductance_grid_precision_level_planar() const {
        return _leakageInductanceGridPrecisionLevelPlanar;
    }
    void Settings::set_leakage_inductance_grid_precision_level_planar(double value) {
        _leakageInductanceGridPrecisionLevelPlanar = value;
    }

    double Settings::get_leakage_inductance_grid_precision_level_wound() const {
        return _leakageInductanceGridPrecisionLevelWound;
    }
    void Settings::set_leakage_inductance_grid_precision_level_wound(double value) {
        _leakageInductanceGridPrecisionLevelWound = value;
    }

    size_t Settings::get_coil_adviser_maximum_number_wires() const {
        return _coilAdviserMaximumNumberWires;
    }
    void Settings::set_coil_adviser_maximum_number_wires(size_t value) {
        _coilAdviserMaximumNumberWires = value;
    }

    bool Settings::get_core_adviser_include_stacks() const {
        return _coreAdviserIncludeStacks;
    }
    void Settings::set_core_adviser_include_stacks(bool value) {
        _coreAdviserIncludeStacks = value;
    }

    bool Settings::get_core_adviser_include_distributed_gaps() const {
        return _coreAdviserIncludeDistributedGaps;
    }
    void Settings::set_core_adviser_include_distributed_gaps(bool value) {
        _coreAdviserIncludeDistributedGaps = value;
    }

    bool Settings::get_core_adviser_include_margin() const {
        return _coreAdviserIncludeMargin;
    }
    void Settings::set_core_adviser_include_margin(bool value) {
        _coreAdviserIncludeMargin = value;
    }

    bool Settings::get_core_adviser_enable_intermediate_pruning() const {
        return _coreAdviserEnableIntermediatePruning;
    }
    void Settings::set_core_adviser_enable_intermediate_pruning(bool value) {
        _coreAdviserEnableIntermediatePruning = value;
    }

    size_t Settings::get_core_adviser_maximum_magnetics_after_filtering() const {
        return _coreAdviserMaximumMagneticsAfterFiltering;
    }
    void Settings::set_core_adviser_maximum_magnetics_after_filtering(size_t value) {
        _coreAdviserMaximumMagneticsAfterFiltering = value;
    }

    bool Settings::get_wire_adviser_include_planar() const {
        return _wireAdviserIncludePlanar;
    }
    void Settings::set_wire_adviser_include_planar(bool value) {
        _wireAdviserIncludePlanar = value;
    }

    bool Settings::get_wire_adviser_include_foil() const {
        return _wireAdviserIncludeFoil;
    }
    void Settings::set_wire_adviser_include_foil(bool value) {
        _wireAdviserIncludeFoil = value;
    }

    bool Settings::get_wire_adviser_include_rectangular() const {
        return _wireAdviserIncludeRectangular;
    }
    void Settings::set_wire_adviser_include_rectangular(bool value) {
        _wireAdviserIncludeRectangular = value;
    }

    bool Settings::get_wire_adviser_include_litz() const {
        return _wireAdviserIncludeLitz;
    }
    void Settings::set_wire_adviser_include_litz(bool value) {
        _wireAdviserIncludeLitz = value;
    }

    bool Settings::get_wire_adviser_include_round() const {
        return _wireAdviserIncludeRound;
    }
    void Settings::set_wire_adviser_include_round(bool value) {
        _wireAdviserIncludeRound = value;
    }

    bool Settings::get_wire_adviser_allow_rectangular_in_toroidal_cores() const {
        return _wireAdviserAllowRectangularInToroidalCores;
    }
    void Settings::set_wire_adviser_allow_rectangular_in_toroidal_cores(bool value) {
        _wireAdviserAllowRectangularInToroidalCores = value;
    }

    bool Settings::get_harmonic_amplitude_threshold_quick_mode() const {
        return _harmonicAmplitudeThresholdQuickMode;
    }
    void Settings::set_harmonic_amplitude_threshold_quick_mode(bool value) {
        _harmonicAmplitudeThresholdQuickMode = value;
    }

    double Settings::get_harmonic_amplitude_threshold() const {
        return _harmonicAmplitudeThreshold;
    }
    void Settings::set_harmonic_amplitude_threshold(double value) {
        _harmonicAmplitudeThreshold = value;
    }

    std::vector<CoreLossesModels> Settings::get_core_losses_model_names() const {
        return _coreLossesModelNames;
    }
    void Settings::set_core_losses_preferred_model_name(CoreLossesModels value) {
        _coreLossesModelNames = {value, CoreLossesModels::PROPRIETARY, CoreLossesModels::LOSS_FACTOR, CoreLossesModels::STEINMETZ, CoreLossesModels::ROSHEN};
    }

    std::string Settings::get_preferred_core_material_ferrite_manufacturer() const {
        return _preferredCoreMaterialFerriteManufacturer;
    }
    void Settings::set_preferred_core_material_ferrite_manufacturer(std::string value) {
        _preferredCoreMaterialFerriteManufacturer = value;
    }

    std::string Settings::get_preferred_core_material_powder_manufacturer() const {
        return _preferredCoreMaterialPowderManufacturer;
    }
    void Settings::set_preferred_core_material_powder_manufacturer(std::string value) {
        _preferredCoreMaterialPowderManufacturer = value;
    }

    bool Settings::get_core_cross_referencer_allow_different_core_material_type() const {
        return _coreCrossReferencerAllowDifferentCoreMaterialType;
    }
    void Settings::set_core_cross_referencer_allow_different_core_material_type(bool value) {
        _coreCrossReferencerAllowDifferentCoreMaterialType = value;
    }

} // namespace OpenMagnetics
