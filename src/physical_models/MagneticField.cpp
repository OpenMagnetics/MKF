#include "physical_models/MagneticField.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "constructive_models/Wire.h"
#include "constructive_models/Magnetic.h"
#include "support/CoilMesher.h"
#include "Models.h"
#include "physical_models/Reluctance.h"
#include "processors/MagneticSimulator.h"
#include "physical_models/InitialPermeability.h"
#include <map>
#include "physical_models/WindingOhmicLosses.h"
#include "support/Settings.h"
#include "support/Utils.h"

#include <algorithm>
#include <cmath>
#include <magic_enum.hpp>

// Some platforms (Emscripten, Apple libc++) don't have std::comp_ellint_1/2
// Use custom implementations from Utils.h via namespace injection
#if defined(__EMSCRIPTEN__) || defined(__APPLE__)
namespace std {
    using OpenMagnetics::comp_ellint_1;
    using OpenMagnetics::comp_ellint_2;
}
#endif
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <streambuf>
#include <vector>
#include <list>
#include <mutex>
#include "support/Exceptions.h"

namespace OpenMagnetics {

namespace {

// Everything the turn-to-point kernel sum reads besides the currents: the mesh (every
// field-point attribute but its value), the winding of each turn filament, the field model,
// the wires and the core geometry the exclusions test. Compared exactly; the hash only
// narrows the search.
struct TurnSumsKey {
    std::vector<double> numbers;
    std::vector<std::string> strings;
    bool operator==(const TurnSumsKey&) const = default;
};

size_t hash_turn_sums_key(const TurnSumsKey& key) {
    // FNV-1a over the bytes.
    uint64_t hash = 1469598103934665603ULL;
    auto mix = [&hash](const void* data, size_t size) {
        auto bytes = static_cast<const unsigned char*>(data);
        for (size_t index = 0; index < size; ++index) {
            hash ^= bytes[index];
            hash *= 1099511628211ULL;
        }
    };
    mix(key.numbers.data(), key.numbers.size() * sizeof(double));
    for (const auto& text : key.strings) {
        mix(text.data(), text.size());
        mix("\0", 1);
    }
    return static_cast<size_t>(hash);
}

struct TurnSumsEntry {
    size_t hash;
    TurnSumsKey key;
    // Per winding, the filament position its currents are measured against; none when the
    // winding carried no current when the sums were built.
    std::vector<std::optional<size_t>> anchorPerWinding;
    // Per filament position: its current over its winding's anchor current.
    std::vector<double> relativeCurrentPerPosition;
    // Per induced point, per winding: (Hx, Hy) of the winding's turns with unit anchor current.
    std::vector<std::vector<std::pair<double, double>>> unitFieldPerInducedPointPerWinding;
    size_t bytes;
};

struct TurnSumsCache {
    std::mutex mutex;
    std::list<TurnSumsEntry> entries;  // most recently used first
    size_t bytes = 0;
    size_t hits = 0;
    size_t misses = 0;
};

TurnSumsCache& turn_sums_cache() {
    static TurnSumsCache cache;
    return cache;
}

size_t turn_sums_entry_bytes(const TurnSumsEntry& entry) {
    size_t bytes = sizeof(TurnSumsEntry) + entry.key.numbers.size() * sizeof(double) +
                   entry.anchorPerWinding.size() * sizeof(std::optional<size_t>) + entry.relativeCurrentPerPosition.size() * sizeof(double);
    for (const auto& text : entry.key.strings) {
        bytes += sizeof(std::string) + text.size();
    }
    for (const auto& perWinding : entry.unitFieldPerInducedPointPerWinding) {
        bytes += sizeof(perWinding) + perWinding.size() * sizeof(std::pair<double, double>);
    }
    return bytes;
}

// Two currents that should be equal up to rounding.
bool currents_agree(double actual, double expected) {
    return std::abs(actual - expected) <= 1e-9 * std::max(std::abs(actual), std::abs(expected));
}

} // namespace

MagneticField::TurnSumsCacheStatistics MagneticField::get_turn_sums_cache_statistics() {
    auto& cache = turn_sums_cache();
    std::lock_guard<std::mutex> lock(cache.mutex);
    return {cache.hits, cache.misses, cache.entries.size(), cache.bytes};
}

void MagneticField::clear_turn_sums_cache() {
    auto& cache = turn_sums_cache();
    std::lock_guard<std::mutex> lock(cache.mutex);
    cache.entries.clear();
    cache.bytes = 0;
    cache.hits = 0;
    cache.misses = 0;
}

SignalDescriptor MagneticField::calculate_magnetic_flux(SignalDescriptor magnetizingCurrent,
                                                          double reluctance,
                                                          double numberTurns) {
    SignalDescriptor magneticFlux;
    Waveform magneticFluxWaveform;
    std::vector<double> magneticFluxData;
    auto compressedMagnetizingCurrentWaveform = magnetizingCurrent.get_waveform().value();

    if (Inputs::is_waveform_sampled(compressedMagnetizingCurrentWaveform)) {
        compressedMagnetizingCurrentWaveform = Inputs::compress_waveform(compressedMagnetizingCurrentWaveform);
    }

    for (auto& datum : compressedMagnetizingCurrentWaveform.get_data()) {
        magneticFluxData.push_back(datum * numberTurns / reluctance);
    }

    if (compressedMagnetizingCurrentWaveform.get_time()) {
        magneticFluxWaveform.set_time(compressedMagnetizingCurrentWaveform.get_time());
    }

    magneticFluxWaveform.set_data(magneticFluxData);
    magneticFlux.set_waveform(magneticFluxWaveform);
    if (magnetizingCurrent.get_harmonics()) {
        auto harmonics = magnetizingCurrent.get_harmonics().value();
        for (size_t harmonicIndex = 0; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
            harmonics.get_mutable_amplitudes()[harmonicIndex] *= numberTurns / reluctance;
        }
        magneticFlux.set_harmonics(harmonics);
    }
    return magneticFlux;
}

SignalDescriptor MagneticField::calculate_magnetic_flux_density(SignalDescriptor magneticFlux,
                                                                  double area) {
    if (!(area > 0)) {
        // B = flux / area: a zero or missing area turns every sample into inf/NaN,
        // which surfaced far away as "Waveform data contains NaN".
        throw std::invalid_argument("calculate_magnetic_flux_density: effective area must be positive, got " +
                                    std::to_string(area));
    }
    SignalDescriptor magneticFluxDensity;
    Waveform magneticFluxDensityWaveform;
    std::vector<double> magneticFluxDensityData;
    auto magneticFluxWaveform = magneticFlux.get_waveform().value();

    if (magneticFluxWaveform.get_time()) {
        magneticFluxDensityWaveform.set_time(magneticFluxWaveform.get_time());
    }

    for (auto& datum : magneticFluxWaveform.get_data()) {
        magneticFluxDensityData.push_back(datum / area);
    }

    magneticFluxDensityWaveform.set_data(magneticFluxDensityData);
    magneticFluxDensity.set_waveform(magneticFluxDensityWaveform);
    if (magneticFlux.get_harmonics()) {
        auto harmonics = magneticFlux.get_harmonics().value();
        for (size_t harmonicIndex = 0; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
            harmonics.get_mutable_amplitudes()[harmonicIndex] /= area;
        }
        magneticFluxDensity.set_harmonics(harmonics);
    }
    magneticFluxDensity.set_processed(
        Inputs::calculate_basic_processed_data(magneticFluxDensityWaveform));

    return magneticFluxDensity;
}

SignalDescriptor MagneticField::calculate_magnetic_field_strength(SignalDescriptor magneticFluxDensity,
                                                                    double initialPermeability) {
    SignalDescriptor magneticFieldStrength;
    Waveform magneticFieldStrengthWaveform;
    std::vector<double> magneticFieldStrengthData;
    auto constants = Constants();
    auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value();

    if (magneticFluxDensityWaveform.get_time()) {
        magneticFieldStrengthWaveform.set_time(magneticFluxDensityWaveform.get_time());
    }

    for (auto& datum : magneticFluxDensityWaveform.get_data()) {
        magneticFieldStrengthData.push_back(datum / (initialPermeability * constants.vacuumPermeability));
    }

    magneticFieldStrengthWaveform.set_data(magneticFieldStrengthData);
    magneticFieldStrength.set_waveform(magneticFieldStrengthWaveform);
    if (magneticFluxDensity.get_harmonics()) {
        auto harmonics = magneticFluxDensity.get_harmonics().value();
        for (size_t harmonicIndex = 0; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
            harmonics.get_mutable_amplitudes()[harmonicIndex] /= (initialPermeability * constants.vacuumPermeability);
        }
        magneticFieldStrength.set_harmonics(harmonics);
    }
    magneticFieldStrength.set_processed(
        Inputs::calculate_basic_processed_data(magneticFieldStrengthWaveform));

    return magneticFieldStrength;
}

std::shared_ptr<MagneticFieldStrengthModel> MagneticField::factory(MagneticFieldStrengthModels modelName) {
    if (modelName == MagneticFieldStrengthModels::BINNS_LAWRENSON) {
        return std::make_shared<MagneticFieldStrengthBinnsLawrensonModel>();
    }
    else if (modelName == MagneticFieldStrengthModels::LAMMERANER) {
        return std::make_shared<MagneticFieldStrengthLammeranerModel>();
    }
    else if (modelName == MagneticFieldStrengthModels::WANG) {
        return std::make_shared<MagneticFieldStrengthWangModel>();
    }
    else if (modelName == MagneticFieldStrengthModels::ALBACH) {
        return std::make_shared<MagneticFieldStrengthAlbach2DModel>();
    }
    else if (modelName == MagneticFieldStrengthModels::DOWELL) {
        return std::make_shared<MagneticFieldStrengthDowellModel>();
    }
    else if (modelName == MagneticFieldStrengthModels::IMAGED_MMF_SHEETS) {
        return std::make_shared<MagneticFieldStrengthImagedMmfSheetsModel>();
    }
    else
        throw ModelNotAvailableException("Unknown Magnetic Field Strength model, available options are: {BINNS_LAWRENSON, LAMMERANER, WANG, ALBACH, DOWELL, IMAGED_MMF_SHEETS}");
}

std::shared_ptr<MagneticFieldStrengthFringingEffectModel> MagneticField::factory(MagneticFieldStrengthFringingEffectModels modelName) {
    if (modelName == MagneticFieldStrengthFringingEffectModels::ALBACH) {
        return std::make_shared<MagneticFieldStrengthAlbachModel>();
    }
    else if (modelName == MagneticFieldStrengthFringingEffectModels::ROSHEN) {
        return std::make_shared<MagneticFieldStrengthRoshenModel>();
    }
    else if (modelName == MagneticFieldStrengthFringingEffectModels::SULLIVAN) {
        return std::make_shared<MagneticFieldStrengthSullivanModel>();
    }
    else
        throw ModelNotAvailableException("Unknown Magnetic Field Strength Fringing Effect model, available options are: {ALBACH, ROSHEN, SULLIVAN}");
}

std::shared_ptr<MagneticFieldStrengthModel> MagneticField::factory() {
    auto defaults = Defaults();
    return factory(defaults.magneticFieldStrengthModelDefault);
}

bool is_inside_inducing_turns(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, Wire* inducingWire) {
    double distanceX = fabs(inducingFieldPoint.get_point()[0] - inducedFieldPoint.get_point()[0]);
    double distanceY = fabs(inducingFieldPoint.get_point()[1] - inducedFieldPoint.get_point()[1]);
    if (inducingWire->get_type() == WireType::ROUND || inducingWire->get_type() == WireType::LITZ) {
        if (hypot(distanceX, distanceY) < inducingWire->get_maximum_outer_width() / 2) {
            return true;
        }
    }
    else {
        if (distanceX < inducingWire->get_maximum_outer_width() / 2 && distanceY < inducingWire->get_maximum_outer_height() / 2) {
            return true;
        }
    }

    return false;
}


bool is_inside_turns(std::vector<Turn> turns, FieldPoint inducedFieldPoint, std::vector<Wire> wires, Magnetic magnetic) {
    for (auto turn : turns) {
        auto windingIndex = magnetic.get_mutable_coil().get_winding_index_by_name(turn.get_winding());
        double distanceX = fabs(turn.get_coordinates()[0] - inducedFieldPoint.get_point()[0]);
        double distanceY = fabs(turn.get_coordinates()[1] - inducedFieldPoint.get_point()[1]);
        if (wires[windingIndex].get_type() == WireType::ROUND || wires[windingIndex].get_type() == WireType::LITZ) {
            if (hypot(distanceX, distanceY) < wires[windingIndex].get_maximum_outer_width() / 2) {
                return true;
            }
        }
        else {
            if (distanceX < wires[windingIndex].get_maximum_outer_width() / 2 && distanceY < wires[windingIndex].get_maximum_outer_height() / 2) {
                return true;
            }
        }
    }

    return false;
}

bool is_inside_core(const FieldPoint& inducedFieldPoint, double coreColumnWidth, double coreWidth, CoreShapeFamily coreShapeFamily) {
    if (coreShapeFamily != CoreShapeFamily::T) {
        return false;
    }
    double radius = sqrt(pow(inducedFieldPoint.get_point()[0], 2) + pow(inducedFieldPoint.get_point()[1], 2));

    if (radius * 1.05 > coreWidth / 2) {
        return false;
    }
    if (radius * 0.95 < (coreWidth / 2 - coreColumnWidth)) {
        return false;
    }
    return true;
}

// Field strength in the gap per ampere-turn of MMF at this frequency: B = MMF / (R_core A_e),
// H_gap = B / mu_0, with the core reluctance at the material's initial permeability at that
// frequency. Used for a single winding, whose current IS the magnetizing current (ABT #1463).
double get_magnetic_field_strength_gap_per_ampere_turn(Magnetic& magnetic, double frequency) {
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory();
    OpenMagnetics::InitialPermeability initial_permeability;
    double initialPermeability = initial_permeability.get_initial_permeability(magnetic.get_mutable_core().resolve_material(), std::nullopt, std::nullopt, frequency);
    double reluctance = reluctanceModel->get_core_reluctance(magnetic.get_core(), initialPermeability).get_core_reluctance();
    double effectiveArea = magnetic.get_core().get_processed_description()->get_effective_parameters().get_effective_area();
    if (!(reluctance > 0) || !(effectiveArea > 0)) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Gap field: core reluctance " + std::to_string(reluctance) +
                                   " H^-1 and effective area " + std::to_string(effectiveArea) + " m2 must both be positive");
    }
    return 1.0 / (reluctance * effectiveArea * Constants().vacuumPermeability);
}

WindingWindowMagneticStrengthFieldPhasorOutput MagneticField::calculate_magnetic_field_strength_field(OperatingPoint operatingPoint, Magnetic magnetic, std::optional<Field> externalInducedField, std::optional<std::vector<int8_t>> customCurrentDirectionPerWinding, std::optional<CoilMesherModels> coilMesherModel) {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    auto includeFringing = settings.get_magnetic_field_include_fringing();

    // ABT #835: DOWELL is a 1-D MMF STAIRCASE, not a 2-D kernel. Its step contribution
    // (fieldStep = I / windingBreadth, applied whole/half/not-at-all by comparing the
    // inducing and induced x) carries NO distance dependence at all — so every filament it
    // is handed counts as a whole extra conductor in the staircase, wherever it sits.
    //
    // The inducing mesh, however, expands each turn into a (2M+1)x(2N+1) lattice of
    // method-of-images filaments, because the 2-D kernels (Lammeraner, Binns-Lawrenson,
    // Albach) NEED those images to satisfy the high-permeability wall boundary condition.
    // Feeding them to Dowell is doubly wrong: the axial images sit at the SAME x as their
    // parent and simply triple every real step, while the inner radial images sit to the
    // left of every real turn and each add a full step — a constant pedestal that dominated
    // the result. On the 12-turn P 3.3/2.6 fixture that read H ~16x high at the first layer
    // and proximity loss 2.27 W against ALBACH/LAMMERANER/BINNS's 0.038 W (~60x, since loss
    // goes as H^2).
    //
    // Dowell's 1-D formulation already embodies the wall boundary condition, so it must see
    // the REAL conductors only. Suppressing mirroring for this model alone brings it to
    // 0.0362 W — within 6% of the 2-D models — while leaving them their images (they drop to
    // 0.0297 W without them, so this cannot be a global change).
    std::optional<SettingsGuard<int>> dowellMirroringGuard;
    if (_magneticFieldStrengthModel == MagneticFieldStrengthModels::DOWELL) {
        dowellMirroringGuard.emplace(settings, &Settings::get_magnetic_field_mirroring_dimension,
                                     &Settings::set_magnetic_field_mirroring_dimension, 0);
    }

    CoilMesher coilMesher; 
    std::vector<Field> inducingFields;
    auto core = magnetic.get_core();

    // A core supplied as functionalDescription-only (e.g. a MAS file saved
    // without processedDescription) has no columns/winding-window geometry yet.
    // The field calculation (and the CoilMesher it invokes via `magnetic`) needs
    // that geometry, so process it here and write it back into `magnetic` —
    // otherwise get_columns() below, and the mesher, would throw
    // CoreNotProcessedException. `magnetic` is taken by value, so this is local.
    if (!core.get_processed_description()) {
        core.process_data();
    }
    if (!core.is_gap_processed()) {
        core.process_gap_or_throw();
    }
    magnetic.set_core(core);
    auto gapping = core.get_functional_description().get_gapping();
    double coreColumnWidth = core.get_columns()[0].get_width();
    auto processedDescription = core.get_processed_description().value();
    double coreWidth = processedDescription.get_width();
    auto coreShapeFamily = core.get_shape_family();

    // ALBACH model uses axisymmetric Biot-Savart (circular current loops), which produces
    // radially symmetric fields. This is incorrect for toroidal cores where turns are placed
    // at specific angular positions -- the field should be localized near the wire, not uniform
    // around the toroid. Fall back to BINNS_LAWRENSON which works in 2D Cartesian and
    // naturally handles the actual (x,y) positions of wire segments.
    if (coreShapeFamily == CoreShapeFamily::T && _magneticFieldStrengthModel == MagneticFieldStrengthModels::ALBACH) {
        _magneticFieldStrengthModel = MagneticFieldStrengthModels::BINNS_LAWRENSON;
        _model = factory(_magneticFieldStrengthModel);
    }

    // Multi-column winding: ALBACH folds r = |x|, modelling every turn as a circular
    // loop around the MAIN column — the wrong topology for turns wound on another
    // column (a lateral leg is not a solid of revolution about the central axis).
    // Fall back to BINNS_LAWRENSON, whose 2D Cartesian filaments plus the per-window
    // mirror images handle any placement.
    bool multiWindowCore = coreShapeFamily != CoreShapeFamily::T && core.get_winding_windows().size() > 1;
    if (multiWindowCore && _magneticFieldStrengthModel == MagneticFieldStrengthModels::ALBACH) {
        _magneticFieldStrengthModel = MagneticFieldStrengthModels::BINNS_LAWRENSON;
        _model = factory(_magneticFieldStrengthModel);
    }

    // MAS excitation convention (2026-09-24): direction by isolation side, not by index.
    std::vector<int8_t> currentDirectionPerWinding;
    if (!customCurrentDirectionPerWinding) {
        currentDirectionPerWinding = CoilMesher::calculate_current_direction_per_winding(magnetic.get_coil());
    }
    else {
        currentDirectionPerWinding = customCurrentDirectionPerWinding.value();
    }

    // Phase of each winding's current per inducing harmonic, aligned with inducingFields.
    std::vector<std::vector<double>> currentPhasePerHarmonicPerWinding;
    std::vector<std::optional<double>> gaugePhasePerHarmonic;
    if (externalInducedField){
        auto aux = coilMesher.generate_mesh_inducing_coil_phasors(magnetic, operatingPoint, settings.get_harmonic_amplitude_threshold(), currentDirectionPerWinding, coilMesherModel);
        // We only process the harmonic that comes from the external field
        for (size_t auxIndex = 0; auxIndex < aux.fieldPerHarmonic.size(); ++auxIndex) {
            if (aux.fieldPerHarmonic[auxIndex].get_frequency() == externalInducedField.value().get_frequency()) {
                inducingFields.push_back(aux.fieldPerHarmonic[auxIndex]);
                currentPhasePerHarmonicPerWinding.push_back(aux.currentPhasePerHarmonicPerWinding[auxIndex]);
                gaugePhasePerHarmonic.push_back(aux.gaugePhasePerHarmonic[auxIndex]);
                break;
            }
        }
    }
    else {
        auto aux = coilMesher.generate_mesh_inducing_coil_phasors(magnetic, operatingPoint, settings.get_harmonic_amplitude_threshold(), currentDirectionPerWinding);
        inducingFields = aux.fieldPerHarmonic;
        currentPhasePerHarmonicPerWinding = aux.currentPhasePerHarmonicPerWinding;
        gaugePhasePerHarmonic = aux.gaugePhasePerHarmonic;
    }
    if (currentPhasePerHarmonicPerWinding.size() != inducingFields.size() || gaugePhasePerHarmonic.size() != inducingFields.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Magnetic field: " + std::to_string(inducingFields.size()) + " inducing harmonics but " +
                                   std::to_string(currentPhasePerHarmonicPerWinding.size()) + " sets of winding current phases");
    }

    if (!magnetic.get_coil().get_turns_description()) {
        throw CoilNotProcessedException("Missing turns description in coil");
    }
    _wirePerWinding = magnetic.get_mutable_coil().get_wires();
    _model->_wirePerWinding = _wirePerWinding;
    
    // Precompute wire dimensions for performance (avoid repeated calculations)
    _wireMaxOuterWidth.clear();
    _wireMaxOuterHeight.clear();
    _wireMaxOuterWidth.reserve(_wirePerWinding.size());
    _wireMaxOuterHeight.reserve(_wirePerWinding.size());
    for (auto& wire : _wirePerWinding) {
        _wireMaxOuterWidth.push_back(wire.get_maximum_outer_width());
        _wireMaxOuterHeight.push_back(wire.get_maximum_outer_height());
    }
    _model->_wireMaxOuterWidth = _wireMaxOuterWidth;
    _model->_wireMaxOuterHeight = _wireMaxOuterHeight;

    // ABT #376: Dowell needs the winding breadth b — the window dimension the layers run along,
    // which for MKF's concentric windows is the window HEIGHT (layers stack radially, the field
    // runs axially). Supplied to every model; only Dowell reads it.
    {
        auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
        if (bobbin.get_processed_description()) {
            auto bobbinWindingWindows = bobbin.get_processed_description()->get_winding_windows();
            if (!bobbinWindingWindows.empty() && bobbinWindingWindows[0].get_height()) {
                _model->_windingWindowBreadth = bobbinWindingWindows[0].get_height().value();
            }
        }
    }
    
    auto turns = magnetic.get_coil().get_turns_description().value();
    std::vector<size_t> windingIndexPerTurn(turns.size());
    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        windingIndexPerTurn[turnIndex] = magnetic.get_mutable_coil().get_winding_index_by_name(turns[turnIndex].get_winding());
    }

    // Set up IMAGED_MMF_SHEETS: turn rectangles/filaments and gap sheets in the image lattice.
    if (_magneticFieldStrengthModel == MagneticFieldStrengthModels::IMAGED_MMF_SHEETS) {
        auto imagedModel = std::dynamic_pointer_cast<MagneticFieldStrengthImagedMmfSheetsModel>(_model);
        if (!imagedModel) {
            throw CalculationException(ErrorCode::CALCULATION_ERROR, "IMAGED_MMF_SHEETS selected but the model instance is of another type");
        }
        auto ohmicLosses = WindingOhmicLosses::calculate_ohmic_losses(magnetic.get_coil(), operatingPoint, defaults.ambientTemperature);
        imagedModel->setup(magnetic, _wirePerWinding, ohmicLosses.get_current_divider_per_turn().value(), currentDirectionPerWinding,
                           operatingPoint.get_excitations_per_winding()[0].get_frequency());
    }

    // Set up ALBACH model if being used
    if (_magneticFieldStrengthModel == MagneticFieldStrengthModels::ALBACH) {
        auto albach2DModel = std::dynamic_pointer_cast<MagneticFieldStrengthAlbach2DModel>(_model);
        if (albach2DModel) {
            albach2DModel->setupFromMagnetic(magnetic, _wirePerWinding);
        }
    }

    std::vector<ComplexField> complexFieldPerHarmonic;
    std::vector<ComplexField> quadratureComplexFieldPerHarmonic;
    for (auto& fieldPerHarmonic : inducingFields) {
        ComplexField field;
        field.set_frequency(fieldPerHarmonic.get_frequency());
        complexFieldPerHarmonic.push_back(field);
        quadratureComplexFieldPerHarmonic.push_back(field);
    }

    std::vector<Field> inducedFields;
    if (externalInducedField) {
        inducedFields.push_back(externalInducedField.value());
    }
    else {
        inducedFields = coilMesher.generate_mesh_induced_coil(magnetic, operatingPoint, settings.get_harmonic_amplitude_threshold());
    }

    // The width-resolved "widthsample" points drive the perpendicular-field integral
    // in the Wang proximity model. That kernel is FEM-validated for the gap-fringing
    // regime (a functional gap dominating the field at a flat conductor); for
    // ungapped cores the inter-conductor perpendicular field it integrates lacks the
    // stack self-screening a real conductor stack exhibits and over-predicts, so the
    // lumped legacy path is kept there: strip the samples when no functional gap exists.
    bool hasFunctionalGap = false;
    for (auto& gap : gapping) {
        if (gap.get_type() == GapType::SUBTRACTIVE || gap.get_type() == GapType::ADDITIVE) {
            hasFunctionalGap = true;
            break;
        }
    }
    if (!hasFunctionalGap) {
        for (auto& inducedField : inducedFields) {
            auto& inducedData = inducedField.get_mutable_data();
            inducedData.erase(std::remove_if(inducedData.begin(), inducedData.end(),
                [](const FieldPoint& point) { return point.get_label() && point.get_label().value() == "widthsample"; }),
                inducedData.end());
        }
    }

    // For ALBACH model, we only compute the turns contribution (air coil field).
    // The gap fringing field is handled by whichever fringing effect model is configured
    // (ALBACH or ROSHEN fringing). This allows users to choose their preferred fringing model.
    // ALBACH supports any number of gaps since the fringing is computed separately.
    bool isAlbach = (_magneticFieldStrengthModel == MagneticFieldStrengthModels::ALBACH);
    
    // Use the configured fringing model (default is ROSHEN, can be overridden to ALBACH)
    std::shared_ptr<MagneticFieldStrengthFringingEffectModel> fringingModel = _fringingEffectModel;

    // Gaps too large for Albach's fitted validity range (xi = gapLength / columnDiameter
    // > ~0.2755), routed through the Roshen conformal model per induced point instead.
    std::vector<CoreGap> albachOutOfRangeGaps;
    MagneticFieldStrengthRoshenModel albachFallbackRoshenModel;

    if (_magneticFieldStrengthFringingEffectModel == MagneticFieldStrengthFringingEffectModels::ALBACH) {
        if (includeFringing) {
            for (auto& gap : gapping) {
                if (gap.get_coordinates().value()[0] < 0) {
                    continue;
                }
                // ABT #832: only functional (SUBTRACTIVE/ADDITIVE) gaps fringe.
                // A RESIDUAL gap is a ground mating surface a few um long; its
                // conformal near-field sampled at a surface point mm away is a
                // modelling artifact, not physics (same doctrine as the
                // width-sample gate below). With residual-gap fringing included
                // a 12-turn P-core read R_ac/R_dc 2.63 at 1 MHz where OMFEM
                // gives 1.263 -- excluding it lands at 1.25.
                if (gap.get_type() != GapType::SUBTRACTIVE && gap.get_type() != GapType::ADDITIVE) {
                    continue;
                }
                // ABT #832 (FEM-arbitrated 2026-08-20): the equivalent-current
                // construction is NOT used, for two independent reasons.
                //  1. Albach's fitted current polynomial (Abb. 9.5) is stated
                //     accurate only for xi = lg/(2rc) < 0.2 and its denominator
                //     has a pole at xi ~ 0.2755; the old validity gate
                //     (denominator > 0) admitted xi up to the pole, where the
                //     equivalent current diverges (ETD24, 2 mm gap, xi = 0.253:
                //     I_eq = 60x the gap MMF -> R_ac/R_dc 53 vs FEM 1.44 at
                //     100 kHz).
                //  2. Even inside the fit's validity the construction needs the
                //     book's full boundary-value treatment (Sect. 9.1.3: the
                //     loop's images in the core) to mean anything; feeding the
                //     equivalent point as a bare wire carrying I = H_g*lg/0.25
                //     (= 4x the gap MMF at small xi) over-predicts fringing
                //     proximity loss 6-9x on an IN-validity gap (ETD24,
                //     0.5 mm gap, xi = 0.063, vs OMFEM).
                // Until the faithful axisymmetric treatment exists, every gap is
                // routed through the Roshen conformal per-point model below --
                // the same path the out-of-validity gaps already took -- which
                // matches OMFEM within 16-35% on the same geometry. Every harmonic
                // gets this gap's field, sized by its own magnetizing flux (ABT #1463).
                albachOutOfRangeGaps.push_back(gap);
            }
        }
    }

    // The field of the turns, summed once for every harmonic. The mesher places each turn's
    // filaments identically at every harmonic and gives filament j of winding w the current
    // base_j * a_w(h): its mesh weight, its turn's current divider and its winding's direction,
    // times the winding's amplitude at that harmonic. Every point-pair kernel is linear in that
    // current and takes no frequency, so winding w's field at an induced point is its anchor
    // filament's current at harmonic h times its field with unit anchor current, and the
    // induced x inducing kernel sum -- the whole cost for a winding of many turns -- runs once
    // instead of once per harmonic. Both properties are checked on the mesh, not assumed: a mesh that breaks either
    // throws. Built on first use: the ALBACH model sums its own way and never needs it. The sums
    // per unit current depend only on the mesh, so they are kept between calls (bounded by
    // Settings::magnetic_field_turn_sums_cache_bytes): re-evaluating a part skips the kernel sum.
    struct TurnFieldSums {
        // Per induced point (in inducedFields order), per winding: (Hx, Hy) of that winding's
        // turns with its anchor filament carrying unit current.
        std::vector<std::vector<std::pair<double, double>>> fieldPerInducedPointPerWinding;
        // Per harmonic, per winding: the anchor filament's current at that harmonic.
        std::vector<std::vector<double>> scalePerHarmonicPerWinding;
    };
    std::optional<TurnFieldSums> turnFieldSums;
    auto build_turn_field_sums = [&]() -> TurnFieldSums {
        size_t numberHarmonics = inducingFields.size();
        size_t numberWindings = magnetic.get_coil().get_functional_description().size();

        // The turn filaments of each harmonic, which must be the same filaments in the same order.
        std::vector<std::vector<size_t>> turnPointsPerHarmonic(numberHarmonics);
        for (size_t harmonicIndex = 0; harmonicIndex < numberHarmonics; ++harmonicIndex) {
            const auto& data = inducingFields[harmonicIndex].get_data();
            for (size_t pointIndex = 0; pointIndex < data.size(); ++pointIndex) {
                if (data[pointIndex].get_turn_index()) {
                    turnPointsPerHarmonic[harmonicIndex].push_back(pointIndex);
                }
            }
        }
        const auto& firstData = inducingFields[0].get_data();
        const auto& firstTurnPoints = turnPointsPerHarmonic[0];
        for (size_t harmonicIndex = 1; harmonicIndex < numberHarmonics; ++harmonicIndex) {
            const auto& data = inducingFields[harmonicIndex].get_data();
            const auto& turnPoints = turnPointsPerHarmonic[harmonicIndex];
            if (turnPoints.size() != firstTurnPoints.size()) {
                throw CalculationException(ErrorCode::CALCULATION_ERROR, "Magnetic field: the coil mesh has " + std::to_string(turnPoints.size()) + " turn filaments at " +
                                           std::to_string(inducingFields[harmonicIndex].get_frequency()) + " Hz but " + std::to_string(firstTurnPoints.size()) + " at the first harmonic");
            }
            for (size_t position = 0; position < turnPoints.size(); ++position) {
                const auto& point = data[turnPoints[position]];
                const auto& firstPoint = firstData[firstTurnPoints[position]];
                if (point.get_point() != firstPoint.get_point() || point.get_turn_index() != firstPoint.get_turn_index()) {
                    throw CalculationException(ErrorCode::CALCULATION_ERROR, "Magnetic field: turn filament " + std::to_string(position) + " of the coil mesh moves between harmonics");
                }
            }
        }
        size_t numberPositions = firstTurnPoints.size();
        std::vector<size_t> windingPerPosition(numberPositions);
        for (size_t position = 0; position < numberPositions; ++position) {
            windingPerPosition[position] = windingIndexPerTurn[static_cast<size_t>(firstData[firstTurnPoints[position]].get_turn_index().value())];
        }
        auto value_at = [&](size_t harmonicIndex, size_t position) {
            return inducingFields[harmonicIndex].get_data()[turnPointsPerHarmonic[harmonicIndex][position]].get_value();
        };

        // The induced points must also be the same at every harmonic.
        const auto& inducedData = inducedFields[0].get_data();
        for (size_t harmonicIndex = 1; harmonicIndex < inducedFields.size(); ++harmonicIndex) {
            const auto& data = inducedFields[harmonicIndex].get_data();
            if (data.size() != inducedData.size()) {
                throw CalculationException(ErrorCode::CALCULATION_ERROR, "Magnetic field: the induced mesh has a different number of points at harmonic " + std::to_string(harmonicIndex));
            }
            for (size_t pointIndex = 0; pointIndex < data.size(); ++pointIndex) {
                if (data[pointIndex].get_point() != inducedData[pointIndex].get_point() || data[pointIndex].get_turn_index() != inducedData[pointIndex].get_turn_index() ||
                    data[pointIndex].get_label() != inducedData[pointIndex].get_label()) {
                    throw CalculationException(ErrorCode::CALCULATION_ERROR, "Magnetic field: induced point " + std::to_string(pointIndex) + " moves between harmonics");
                }
            }
        }

        // What the sum depends on besides the currents.
        TurnSumsKey key;
        auto add_point = [&key](const FieldPoint& point) {
            const auto& coordinates = point.get_point();
            key.numbers.push_back(static_cast<double>(coordinates.size()));
            key.numbers.insert(key.numbers.end(), coordinates.begin(), coordinates.end());
            for (const auto& optionalNumber : {point.get_turn_index() ? std::optional<double>(static_cast<double>(point.get_turn_index().value())) : std::nullopt,
                                               point.get_turn_length(), point.get_rotation()}) {
                key.numbers.push_back(optionalNumber ? 1 : 0);
                key.numbers.push_back(optionalNumber.value_or(0));
            }
            key.strings.push_back(point.get_label() ? "1" + point.get_label().value() : "0");
        };
        key.numbers = {static_cast<double>(_magneticFieldStrengthModel), multiWindowCore ? 1.0 : 0.0, coreColumnWidth, coreWidth,
                       static_cast<double>(coreShapeFamily), _model->_windingWindowBreadth, static_cast<double>(numberWindings)};
        for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
            json wireJson;
            to_json(wireJson, _wirePerWinding[windingIndex]);
            key.strings.push_back(wireJson.dump());
            key.numbers.push_back(_wireMaxOuterWidth[windingIndex]);
            key.numbers.push_back(_wireMaxOuterHeight[windingIndex]);
        }
        key.numbers.push_back(static_cast<double>(numberPositions));
        for (size_t position = 0; position < numberPositions; ++position) {
            add_point(firstData[firstTurnPoints[position]]);
            key.numbers.push_back(static_cast<double>(windingPerPosition[position]));
        }
        key.numbers.push_back(static_cast<double>(inducedData.size()));
        for (const auto& inducedPoint : inducedData) {
            add_point(inducedPoint);
        }
        size_t keyHash = hash_turn_sums_key(key);

        // The anchor currents of this operating point, per harmonic and winding.
        auto scales_for = [&](const std::vector<std::optional<size_t>>& anchorPerWinding) {
            std::vector<std::vector<double>> scalePerHarmonicPerWinding(numberHarmonics, std::vector<double>(numberWindings, 0.0));
            for (size_t harmonicIndex = 0; harmonicIndex < numberHarmonics; ++harmonicIndex) {
                for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
                    if (anchorPerWinding[windingIndex]) {
                        scalePerHarmonicPerWinding[harmonicIndex][windingIndex] = value_at(harmonicIndex, anchorPerWinding[windingIndex].value());
                    }
                }
            }
            return scalePerHarmonicPerWinding;
        };
        // Do this operating point's currents stand to the anchors as the stored ones did?
        auto currents_fit = [&](const TurnSumsEntry& entry) {
            for (size_t position = 0; position < numberPositions; ++position) {
                const auto& anchor = entry.anchorPerWinding[windingPerPosition[position]];
                for (size_t harmonicIndex = 0; harmonicIndex < numberHarmonics; ++harmonicIndex) {
                    double actual = value_at(harmonicIndex, position);
                    double expected = anchor ? entry.relativeCurrentPerPosition[position] * value_at(harmonicIndex, anchor.value()) : 0.0;
                    if (!currents_agree(actual, expected)) {
                        return false;
                    }
                }
            }
            return true;
        };

        size_t cacheLimit = settings.get_magnetic_field_turn_sums_cache_bytes();
        auto& cache = turn_sums_cache();
        if (cacheLimit > 0) {
            std::lock_guard<std::mutex> lock(cache.mutex);
            for (auto entry = cache.entries.begin(); entry != cache.entries.end(); ++entry) {
                if (entry->hash != keyHash || !(entry->key == key)) {
                    continue;
                }
                if (currents_fit(*entry)) {
                    cache.entries.splice(cache.entries.begin(), cache.entries, entry);
                    ++cache.hits;
                    return {entry->unitFieldPerInducedPointPerWinding, scales_for(entry->anchorPerWinding)};
                }
                // Same mesh, currents distributed differently: summed afresh and replaced below.
                cache.bytes -= entry->bytes;
                cache.entries.erase(entry);
                break;
            }
            ++cache.misses;
        }

        // Per winding: the filament and harmonic carrying its largest current. Every harmonic's
        // currents must be proportional to that harmonic's.
        TurnSumsEntry entry;
        entry.hash = keyHash;
        entry.anchorPerWinding.assign(numberWindings, std::nullopt);
        entry.relativeCurrentPerPosition.assign(numberPositions, 0.0);
        std::vector<size_t> referenceHarmonicPerWinding(numberWindings, 0);
        for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
            double largest = 0;
            for (size_t harmonicIndex = 0; harmonicIndex < numberHarmonics; ++harmonicIndex) {
                for (size_t position = 0; position < numberPositions; ++position) {
                    if (windingPerPosition[position] == windingIndex && std::abs(value_at(harmonicIndex, position)) > largest) {
                        largest = std::abs(value_at(harmonicIndex, position));
                        referenceHarmonicPerWinding[windingIndex] = harmonicIndex;
                        entry.anchorPerWinding[windingIndex] = position;
                    }
                }
            }
        }
        for (size_t position = 0; position < numberPositions; ++position) {
            const auto& anchor = entry.anchorPerWinding[windingPerPosition[position]];
            if (anchor) {
                size_t reference = referenceHarmonicPerWinding[windingPerPosition[position]];
                entry.relativeCurrentPerPosition[position] = value_at(reference, position) / value_at(reference, anchor.value());
            }
        }
        if (!currents_fit(entry)) {
            throw CalculationException(ErrorCode::CALCULATION_ERROR, "Magnetic field: the currents of a winding's turn filaments are not proportional across harmonics, "
                                       "so one field sum cannot serve every harmonic");
        }

        // The kernel sum with unit anchor current, with the same exclusions the per-harmonic
        // loop applies to a turn.
        entry.unitFieldPerInducedPointPerWinding.assign(inducedData.size(), std::vector<std::pair<double, double>>(numberWindings, {0.0, 0.0}));
        for (size_t inducedIndex = 0; inducedIndex < inducedData.size(); ++inducedIndex) {
            const auto& inducedFieldPoint = inducedData[inducedIndex];
            // Width samples receive no turn field (the fringing-only rule of the loop below).
            if (inducedFieldPoint.get_label() && inducedFieldPoint.get_label().value() == "widthsample") {
                continue;
            }
            bool inducedPointInsideCore = is_inside_core(inducedFieldPoint, coreColumnWidth, coreWidth, coreShapeFamily);
            if (!inducedFieldPoint.get_turn_index() && inducedPointInsideCore) {
                continue;
            }
            auto& fieldPerWinding = entry.unitFieldPerInducedPointPerWinding[inducedIndex];
            for (size_t position = 0; position < numberPositions; ++position) {
                size_t windingIndex = windingPerPosition[position];
                if (!entry.anchorPerWinding[windingIndex]) {
                    continue;
                }
                size_t reference = referenceHarmonicPerWinding[windingIndex];
                const auto& inducingFieldPoint = inducingFields[reference].get_data()[turnPointsPerHarmonic[reference][position]];
                if (multiWindowCore && inducingFieldPoint.get_point()[0] * inducedFieldPoint.get_point()[0] < 0) {
                    continue;
                }
                if (inducedFieldPoint.get_turn_index() && inducedFieldPoint.get_turn_index().value() == inducingFieldPoint.get_turn_index().value()) {
                    continue;
                }
                auto [inducedFieldX, inducedFieldY] = _model->get_magnetic_field_strength_components_between_two_points(inducingFieldPoint, inducedFieldPoint, windingIndex);
                if (std::isnan(inducedFieldX) || std::isnan(inducedFieldY)) {
                    throw NaNResultException("NaN found in magnetic field calculation");
                }
                fieldPerWinding[windingIndex].first += inducedFieldX;
                fieldPerWinding[windingIndex].second += inducedFieldY;
            }
        }
        for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
            if (!entry.anchorPerWinding[windingIndex]) {
                continue;
            }
            double anchorCurrent = value_at(referenceHarmonicPerWinding[windingIndex], entry.anchorPerWinding[windingIndex].value());
            for (auto& fieldPerWinding : entry.unitFieldPerInducedPointPerWinding) {
                fieldPerWinding[windingIndex].first /= anchorCurrent;
                fieldPerWinding[windingIndex].second /= anchorCurrent;
            }
        }

        TurnFieldSums sums{entry.unitFieldPerInducedPointPerWinding, scales_for(entry.anchorPerWinding)};
        if (cacheLimit > 0) {
            entry.key = std::move(key);
            entry.bytes = turn_sums_entry_bytes(entry);
            std::lock_guard<std::mutex> lock(cache.mutex);
            // A mesh larger than the whole budget is not kept, rather than emptying the cache for it.
            if (entry.bytes <= cacheLimit) {
                cache.bytes += entry.bytes;
                cache.entries.push_front(std::move(entry));
                while (cache.bytes > cacheLimit) {
                    cache.bytes -= cache.entries.back().bytes;
                    cache.entries.pop_back();
                }
            }
        }
        return sums;
    };

    for (size_t harmonicIndex = 0; harmonicIndex < inducingFields.size(); ++harmonicIndex){
        std::vector<ComplexFieldPoint> fieldPoints;
        std::vector<ComplexFieldPoint> quadratureFieldPoints;

        if (inducedFields[harmonicIndex].get_data().size() == 0) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Empty complexField");
        }

        // Every kernel is LINEAR in the inducing current, so a turn of winding k with complex
        // current value * exp(j phase_k) contributes (Hx, Hy) * cos(phase_k) in phase and
        // (Hx, Hy) * sin(phase_k) in quadrature, from ONE kernel evaluation. The gauge winding
        // has phase exactly 0 (cos 1, sin 0), so a single winding or windings in exact antiphase
        // reproduce the amplitude-only in-phase field bit for bit.
        const auto& currentPhasePerWinding = currentPhasePerHarmonicPerWinding[harmonicIndex];
        std::vector<double> inPhaseFactorPerWinding(currentPhasePerWinding.size());
        std::vector<double> quadratureFactorPerWinding(currentPhasePerWinding.size());
        bool anyQuadratureCurrent = false;
        for (size_t windingIndex = 0; windingIndex < currentPhasePerWinding.size(); ++windingIndex) {
            inPhaseFactorPerWinding[windingIndex] = std::cos(currentPhasePerWinding[windingIndex]);
            quadratureFactorPerWinding[windingIndex] = std::sin(currentPhasePerWinding[windingIndex]);
            if (quadratureFactorPerWinding[windingIndex] != 0) {
                anyQuadratureCurrent = true;
            }
        }

        // The gap fringing field of THIS harmonic (ABT #1463): every harmonic of the magnetizing
        // flux drives its own gap field, on its own phase. The gap carries the core's magnetizing
        // MMF, which by Ampere around the core path is sum_k N_k i_k = R_core * Phi: the load
        // currents cancel there by definition. So the magnetizing MMF the winding currents "already
        // carry" IS their net MMF M_w, and "M_w + N_p I_mag - (the part of M_w that is magnetizing)"
        // is exactly N_p I_mag = R_core Phi: no threshold, no decomposition of the currents. The
        // flux itself comes from the source that is exact for the input:
        //  - one winding: its current is the magnetizing current (nothing can cancel it), so
        //    M_h = c N I_h (MAS amplitude, DFT phase against the gauge -- phase 0 when it is the
        //    gauge) and H_gap,h = M_h / (R_core(f_h) A_e mu_0);
        //  - several windings: Faraday on winding 0, v_0 = N_0 dPhi/dt (dot convention; Phi counts
        //    positive for current into the dot, which is the turn-field frame: filament current =
        //    direction x MAS current = into-dot current). With x(t) = Re(X e^{j w t}),
        //    B_h = V_0,h / (j w_h N_0 A_e) and H_gap,h = B_h / mu_0. For a flyback the windings'
        //    currents carry the magnetizing MMF and Faraday reproduces their net MMF; for a
        //    transformer given ideal (exactly cancelling) currents Faraday still gives the
        //    magnetizing field, which the net MMF (zero) would lose. MKF's derived
        //    magnetizingCurrent is NOT used: for FLYBACK_PRIMARY/UNIPOLAR_TRIANGULAR labels it is a
        //    triangle rebuilt from the current's duty cycle that ignores a DCM dead time (h1 on a
        //    DCM flyback: 17% low and 29 degrees off the windings' net MMF).
        // The DC bin drives no eddy loss and gets no gap field (as before).
        // IMAGED_MMF_SHEETS (below) still sizes and phases its gap sheets from the derived
        // magnetizingCurrent, and so shares that DCM defect; left as is.
        struct GapSource {
            double magnitude;
            double inPhaseFactor;
            double quadratureFactor;
        };
        std::optional<GapSource> gapSource;
        auto get_gap_source = [&]() -> GapSource {
            if (gapSource) {
                return gapSource.value();
            }
            double harmonicFrequency = inducingFields[harmonicIndex].get_frequency();
            if (!(harmonicFrequency > 0)) {
                gapSource = GapSource{0, 0, 0};
                return gapSource.value();
            }
            if (!gaugePhasePerHarmonic[harmonicIndex]) {
                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Gap field at " + std::to_string(harmonicFrequency) +
                                           " Hz: no winding current carries this harmonic, so there is no phase reference");
            }
            const auto& functionalDescription = magnetic.get_coil().get_functional_description();
            const auto& excitations = operatingPoint.get_excitations_per_winding();
            if (excitations.empty() || functionalDescription.empty()) {
                throw InvalidInputException(ErrorCode::MISSING_DATA, "Gap field: no windings or no excitations");
            }
            std::complex<double> fieldPhasor;
            if (functionalDescription.size() == 1) {
                auto current = excitations[0].get_current();
                if (!current || !current->get_harmonics()) {
                    throw InvalidInputException(ErrorCode::MISSING_DATA, "Gap field at " + std::to_string(harmonicFrequency) + " Hz: the winding has no current harmonics");
                }
                const auto harmonics = current->get_harmonics().value();
                double amplitude = 0;
                bool found = false;
                for (size_t index = 0; index < harmonics.get_frequencies().size(); ++index) {
                    if (harmonics.get_frequencies()[index] == harmonicFrequency) {
                        amplitude = harmonics.get_amplitudes()[index];
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    throw CalculationException(ErrorCode::CALCULATION_ERROR, "Gap field at " + std::to_string(harmonicFrequency) +
                                               " Hz: the winding current's harmonics do not list the inducing harmonic");
                }
                if (amplitude == 0) {
                    gapSource = GapSource{0, 0, 0};
                    return gapSource.value();
                }
                double numberTurns = static_cast<double>(functionalDescription[0].get_number_turns());
                double magnetomotiveForce = currentDirectionPerWinding.at(0) * numberTurns * amplitude;
                fieldPhasor = std::polar(magnetomotiveForce * get_magnetic_field_strength_gap_per_ampere_turn(magnetic, harmonicFrequency),
                                         currentPhasePerHarmonicPerWinding[harmonicIndex].at(0));
            }
            else {
                auto voltage = excitations[0].get_voltage();
                if (!voltage || !voltage->get_waveform()) {
                    throw InvalidInputException(ErrorCode::MISSING_DATA, "Gap field at " + std::to_string(harmonicFrequency) +
                                                " Hz: with several windings the magnetizing flux comes from winding 0's voltage (Faraday), and it has no voltage waveform");
                }
                double effectiveArea = magnetic.get_core().get_processed_description()->get_effective_parameters().get_effective_area();
                double numberTurns = static_cast<double>(functionalDescription[0].get_number_turns());
                if (!(effectiveArea > 0) || !(numberTurns > 0)) {
                    throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Gap field: effective area " + std::to_string(effectiveArea) +
                                               " m2 and winding 0 turns " + std::to_string(numberTurns) + " must both be positive");
                }
                auto voltagePhasor = CoilMesher::calculate_harmonic_phasor(voltage->get_waveform().value(), excitations[0].get_frequency(), harmonicFrequency);
                double angularFrequency = 2 * std::numbers::pi * harmonicFrequency;
                auto fluxDensityPhasor = voltagePhasor / (std::complex<double>(0, 1) * angularFrequency * numberTurns * effectiveArea);
                fieldPhasor = fluxDensityPhasor / Constants().vacuumPermeability * std::polar(1.0, -gaugePhasePerHarmonic[harmonicIndex].value());
            }
            double magnitude = std::abs(fieldPhasor);
            if (!(magnitude > 0)) {
                gapSource = GapSource{0, 0, 0};
                return gapSource.value();
            }
            gapSource = GapSource{magnitude, fieldPhasor.real() / magnitude, fieldPhasor.imag() / magnitude};
            return gapSource.value();
        };

        // IMAGED_MMF_SHEETS: the whole window field (turns and gap sheets) at every induced point,
        // width samples included -- the sheets only balance the turns when both are present.
        if (_magneticFieldStrengthModel == MagneticFieldStrengthModels::IMAGED_MMF_SHEETS) {
            auto imagedModel = std::dynamic_pointer_cast<MagneticFieldStrengthImagedMmfSheetsModel>(_model);
            double harmonicFrequency = inducingFields[harmonicIndex].get_frequency();
            size_t numberWindings = magnetic.get_coil().get_functional_description().size();
            // Each winding's current at this harmonic, as the CoilMesher reads it (a winding
            // without this harmonic carries none of it).
            std::vector<double> inPhaseCurrentPerWinding(numberWindings, 0.0);
            std::vector<double> quadratureCurrentPerWinding(numberWindings, 0.0);
            for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
                auto current = operatingPoint.get_excitations_per_winding()[windingIndex].get_current();
                if (!current || !current->get_harmonics()) {
                    throw InvalidInputException(ErrorCode::MISSING_DATA, "IMAGED_MMF_SHEETS: winding " + std::to_string(windingIndex) + " has no current harmonics");
                }
                auto harmonics = current->get_harmonics().value();
                for (size_t k = 0; k < harmonics.get_frequencies().size(); ++k) {
                    if (std::abs(harmonics.get_frequencies()[k] - harmonicFrequency) <= 1e-9 * harmonicFrequency) {
                        inPhaseCurrentPerWinding[windingIndex] = harmonics.get_amplitudes()[k] * inPhaseFactorPerWinding[windingIndex];
                        quadratureCurrentPerWinding[windingIndex] = harmonics.get_amplitudes()[k] * quadratureFactorPerWinding[windingIndex];
                        break;
                    }
                }
            }
            // The gap sheets carry the magnetizing current of this harmonic, on its phase.
            double gapInPhaseCurrent = 0;
            double gapQuadratureCurrent = 0;
            if (includeFringing && !imagedModel->get_gap_sources().empty()) {
                if (!operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current()) {
                    auto magnetizingInductance = MagneticSimulator().calculate_magnetizing_inductance(operatingPoint, magnetic);
                    auto includeDcCurrent = Inputs::include_dc_offset_into_magnetizing_current(operatingPoint, magnetic.get_turns_ratios());
                    auto magnetizingCurrent = Inputs::calculate_magnetizing_current(operatingPoint.get_mutable_excitations_per_winding()[0],
                                                                                   resolve_dimensional_values(magnetizingInductance.get_magnetizing_inductance()),
                                                                                   true, includeDcCurrent,
                                                                                   operatingPoint.get_excitations_per_winding().size() > 1);
                    operatingPoint.get_mutable_excitations_per_winding()[0].set_magnetizing_current(magnetizingCurrent);
                }
                auto magnetizingCurrentSignal = operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current().value();
                if (!magnetizingCurrentSignal.get_waveform()) {
                    throw InvalidInputException(ErrorCode::MISSING_DATA, "IMAGED_MMF_SHEETS: the magnetizing current has no waveform");
                }
                auto magnetizingWaveform = magnetizingCurrentSignal.get_waveform().value();
                auto phasor = CoilMesher::calculate_harmonic_phasor(magnetizingWaveform, operatingPoint.get_excitations_per_winding()[0].get_frequency(), harmonicFrequency);
                double waveformScale = 0;
                for (auto value : magnetizingWaveform.get_data()) {
                    waveformScale = std::max(waveformScale, std::abs(value));
                }
                if (std::abs(phasor) > 1e-9 * waveformScale) {
                    if (!gaugePhasePerHarmonic[harmonicIndex]) {
                        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "IMAGED_MMF_SHEETS at " + std::to_string(harmonicFrequency) +
                                                   " Hz: no winding current carries this harmonic, so there is no phase reference");
                    }
                    double magnetizingPhase = std::arg(phasor) - gaugePhasePerHarmonic[harmonicIndex].value();
                    gapInPhaseCurrent = std::abs(phasor) * std::cos(magnetizingPhase);
                    gapQuadratureCurrent = std::abs(phasor) * std::sin(magnetizingPhase);
                }
            }
            for (const auto& inducedFieldPoint : inducedFields[harmonicIndex].get_data()) {
                double x = inducedFieldPoint.get_point()[0];
                double y = inducedFieldPoint.get_point()[1];
                auto excludedTurn = inducedFieldPoint.get_turn_index();
                auto [turnsX, turnsY] = imagedModel->turns_field(x, y, inPhaseCurrentPerWinding, excludedTurn);
                double quadratureX = 0, quadratureY = 0;
                if (anyQuadratureCurrent) {
                    std::tie(quadratureX, quadratureY) = imagedModel->turns_field(x, y, quadratureCurrentPerWinding, excludedTurn);
                }
                double totalX = turnsX, totalY = turnsY;
                if (gapInPhaseCurrent != 0 || gapQuadratureCurrent != 0) {
                    auto [gapsX, gapsY] = imagedModel->gaps_field(x, y);
                    totalX += gapsX * gapInPhaseCurrent;
                    totalY += gapsY * gapInPhaseCurrent;
                    quadratureX += gapsX * gapQuadratureCurrent;
                    quadratureY += gapsY * gapQuadratureCurrent;
                }
                if (std::isnan(totalX) || std::isnan(totalY) || std::isnan(quadratureX) || std::isnan(quadratureY)) {
                    throw NaNResultException("NaN found in IMAGED_MMF_SHEETS magnetic field calculation");
                }
                ComplexFieldPoint complexFieldPoint;
                complexFieldPoint.set_point(inducedFieldPoint.get_point());
                complexFieldPoint.set_real(totalX);
                complexFieldPoint.set_imaginary(totalY);
                if (excludedTurn) {
                    complexFieldPoint.set_turn_index(excludedTurn.value());
                }
                if (inducedFieldPoint.get_label()) {
                    complexFieldPoint.set_label(inducedFieldPoint.get_label().value());
                }
                fieldPoints.push_back(complexFieldPoint);
                ComplexFieldPoint quadratureFieldPoint(complexFieldPoint);
                quadratureFieldPoint.set_real(quadratureX);
                quadratureFieldPoint.set_imaginary(quadratureY);
                quadratureFieldPoints.push_back(quadratureFieldPoint);
            }
            complexFieldPerHarmonic[harmonicIndex].set_data(fieldPoints);
            quadratureComplexFieldPerHarmonic[harmonicIndex].set_data(quadratureFieldPoints);
            continue;
        }

        // For ALBACH model, use a more efficient approach that calculates
        // the total field from all turns at once for each induced point.
        if (_magneticFieldStrengthModel == MagneticFieldStrengthModels::ALBACH) {
            auto albach2DModel = std::dynamic_pointer_cast<MagneticFieldStrengthAlbach2DModel>(_model);
            if (albach2DModel) {
                // Update turn currents based on harmonic data
                auto& harmonicData = inducingFields[harmonicIndex].get_data();
                std::vector<double> turnCurrents(turns.size(), 0.0);
                std::vector<double> quadratureTurnCurrents(turns.size(), 0.0);
                
                // Collect fringing field inducing points (those without turn_index)
                std::vector<FieldPoint> fringingPoints;
                
                for (auto& inducingPoint : harmonicData) {
                    if (inducingPoint.get_turn_index()) {
                        size_t turnIdx = inducingPoint.get_turn_index().value();
                        if (turnIdx < turnCurrents.size()) {
                            size_t turnWindingIndex = windingIndexPerTurn[turnIdx];
                            turnCurrents[turnIdx] = inducingPoint.get_value() * inPhaseFactorPerWinding[turnWindingIndex];
                            quadratureTurnCurrents[turnIdx] = inducingPoint.get_value() * quadratureFactorPerWinding[turnWindingIndex];
                        }
                    } else {
                        // This is a fringing field equivalent point (from ALBACH fringing model)
                        fringingPoints.push_back(inducingPoint);
                    }
                }
                albach2DModel->updateTurnCurrents(turnCurrents);

                // Update skin depths for frequency-dependent current distribution (Wang 2018)
                // At high frequency, current concentrates at conductor edges
                double frequency = inducingFields[harmonicIndex].get_frequency();
                if (frequency > 0 && !_wirePerWinding.empty()) {
                    // Use the first wire to calculate a representative skin depth
                    // In practice, all wires in a winding should have similar material
                    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(
                        _wirePerWinding[0], frequency, operatingPoint.get_conditions().get_ambient_temperature());
                    albach2DModel->updateSkinDepths(skinDepth);
                }
                
                // Create a model for computing fringing field contribution from equivalent current loops (ALBACH model)
                // Using BINNS_LAWRENSON to compute the field from the equivalent current loops
                auto fringingFieldModel = factory(MagneticFieldStrengthModels::BINNS_LAWRENSON);
                
                // For ROSHEN and SULLIVAN fringing, we need the gap field strength
                // (both use direct gap-to-point calculation via get_magnetic_field_strength_between_gap_and_point).
                // Also needed when ALBACH fringing routed out-of-validity gaps to Roshen.
                double magneticFieldStrengthGap = 0;
                if ((_magneticFieldStrengthFringingEffectModel == MagneticFieldStrengthFringingEffectModels::ROSHEN ||
                     _magneticFieldStrengthFringingEffectModel == MagneticFieldStrengthFringingEffectModels::SULLIVAN ||
                     !albachOutOfRangeGaps.empty()) && includeFringing && hasFunctionalGap) {
                    magneticFieldStrengthGap = get_gap_source().magnitude;
                }
                
                // Calculate field at each induced point directly from all turns
                // Quadrature part of the fringing field per emitted point (the in-phase part is
                // added in place), replayed by the quadrature sweep below.
                std::vector<std::pair<double, double>> fringingQuadraturePerPoint;
                for (auto& inducedFieldPoint : inducedFields[harmonicIndex].get_data()) {
                    // Skip points inside the core
                    if (is_inside_core(inducedFieldPoint, coreColumnWidth, coreWidth, coreShapeFamily)) {
                        continue;
                    }

                    // Width samples feed the Wang perpendicular-field integral, whose
                    // FEM validation covers ONLY the gap-fringing field. Turn fields at
                    // these points are filament-discretization artifacts: a neighbouring
                    // foil/planar sheet reduced to filaments reads kA/m of perpendicular
                    // field where the real co-extensive sheet produces ~none (ABT #182,
                    // 208 W on a 4-turn foil whose Dowell loss is ~2 W), and the induced
                    // turn's own filaments add subdivision noise on top. Inter-conductor
                    // proximity stays with the lumped labeled points.
                    bool fringingOnlyPoint = inducedFieldPoint.get_label() &&
                                             inducedFieldPoint.get_label().value() == "widthsample";
                    ComplexFieldPoint complexFieldPoint;
                    if (fringingOnlyPoint) {
                        complexFieldPoint.set_real(0);
                        complexFieldPoint.set_imaginary(0);
                        complexFieldPoint.set_point(inducedFieldPoint.get_point());
                        if (inducedFieldPoint.get_turn_index()) {
                            complexFieldPoint.set_turn_index(inducedFieldPoint.get_turn_index().value());
                        }
                        complexFieldPoint.set_label(inducedFieldPoint.get_label().value());
                    }
                    else {
                        complexFieldPoint = albach2DModel->calculateTotalFieldAtPoint(inducedFieldPoint);
                    }
                    
                    // Add fringing field contribution based on configured fringing model. The
                    // fringing field carries the magnetizing flux's phase (see
                    // get_gap_source): cos in phase, sin in quadrature.
                    double fringingQuadratureX = 0;
                    double fringingQuadratureY = 0;
                    auto add_fringing = [&](const ComplexFieldPoint& contribution) {
                        double inPhaseFactor = get_gap_source().inPhaseFactor;
                        double quadratureFactor = get_gap_source().quadratureFactor;
                        complexFieldPoint.set_real(complexFieldPoint.get_real() + contribution.get_real() * inPhaseFactor);
                        complexFieldPoint.set_imaginary(complexFieldPoint.get_imaginary() + contribution.get_imaginary() * inPhaseFactor);
                        fringingQuadratureX += contribution.get_real() * quadratureFactor;
                        fringingQuadratureY += contribution.get_imaginary() * quadratureFactor;
                    };
                    if (includeFringing && hasFunctionalGap && get_gap_source().magnitude > 0) {
                        if (_magneticFieldStrengthFringingEffectModel == MagneticFieldStrengthFringingEffectModels::ALBACH) {
                            // ALBACH fringing: use equivalent current loops
                            for (auto& fringingPoint : fringingPoints) {
                                auto fringingContrib = fringingFieldModel->get_magnetic_field_strength_between_two_points(fringingPoint, inducedFieldPoint);
                                add_fringing(fringingContrib);
                            }
                            // Gaps beyond Albach's fitted validity: Roshen conformal model
                            for (auto& gap : albachOutOfRangeGaps) {
                                auto fringingContrib = albachFallbackRoshenModel.get_magnetic_field_strength_between_gap_and_point(gap, magneticFieldStrengthGap, inducedFieldPoint);
                                add_fringing(fringingContrib);
                            }
                        } else if (_magneticFieldStrengthFringingEffectModel == MagneticFieldStrengthFringingEffectModels::ROSHEN ||
                                   _magneticFieldStrengthFringingEffectModel == MagneticFieldStrengthFringingEffectModels::SULLIVAN) {
                            // ROSHEN and SULLIVAN fringing: compute field directly from each gap
                            for (auto& gap : gapping) {
                                if (gap.get_coordinates().value()[0] < 0) {
                                    continue;
                                }
                                // ABT #832: residual (mating-surface) gaps do not contribute
                                // fringing loss -- see the matching gate in the ALBACH branch.
                                if (gap.get_type() != GapType::SUBTRACTIVE && gap.get_type() != GapType::ADDITIVE) {
                                    continue;
                                }
                                auto fringingContrib = _fringingEffectModel->get_magnetic_field_strength_between_gap_and_point(gap, magneticFieldStrengthGap, inducedFieldPoint);
                                add_fringing(fringingContrib);
                            }
                        }
                    }
                    
                    if (std::isnan(complexFieldPoint.get_real()) || std::isnan(complexFieldPoint.get_imaginary())) {
                        throw NaNResultException("NaN found in ALBACH magnetic field calculation");
                    }
                    
                    fieldPoints.push_back(complexFieldPoint);
                    fringingQuadraturePerPoint.push_back({fringingQuadratureX, fringingQuadratureY});
                }

                // Quadrature sweep: the turns' quadrature currents plus the quadrature part of the
                // gap fringing field (the magnetizing current's phase) kept from the sweep above.
                if (anyQuadratureCurrent) {
                    albach2DModel->updateTurnCurrents(quadratureTurnCurrents);
                }
                size_t quadraturePointIndex = 0;
                for (auto& inducedFieldPoint : inducedFields[harmonicIndex].get_data()) {
                    if (is_inside_core(inducedFieldPoint, coreColumnWidth, coreWidth, coreShapeFamily)) {
                        continue;
                    }
                    bool fringingOnlyPoint = inducedFieldPoint.get_label() &&
                                             inducedFieldPoint.get_label().value() == "widthsample";
                    ComplexFieldPoint quadratureFieldPoint;
                    if (fringingOnlyPoint || !anyQuadratureCurrent) {
                        quadratureFieldPoint.set_real(0);
                        quadratureFieldPoint.set_imaginary(0);
                        quadratureFieldPoint.set_point(inducedFieldPoint.get_point());
                        if (inducedFieldPoint.get_turn_index()) {
                            quadratureFieldPoint.set_turn_index(inducedFieldPoint.get_turn_index().value());
                        }
                        if (inducedFieldPoint.get_label()) {
                            quadratureFieldPoint.set_label(inducedFieldPoint.get_label().value());
                        }
                    }
                    else {
                        quadratureFieldPoint = albach2DModel->calculateTotalFieldAtPoint(inducedFieldPoint);
                    }
                    quadratureFieldPoint.set_real(quadratureFieldPoint.get_real() + fringingQuadraturePerPoint[quadraturePointIndex].first);
                    quadratureFieldPoint.set_imaginary(quadratureFieldPoint.get_imaginary() + fringingQuadraturePerPoint[quadraturePointIndex].second);
                    ++quadraturePointIndex;
                    if (std::isnan(quadratureFieldPoint.get_real()) || std::isnan(quadratureFieldPoint.get_imaginary())) {
                        throw NaNResultException("NaN found in ALBACH quadrature magnetic field calculation");
                    }
                    quadratureFieldPoints.push_back(quadratureFieldPoint);
                }
                complexFieldPerHarmonic[harmonicIndex].set_data(fieldPoints);
                quadratureComplexFieldPerHarmonic[harmonicIndex].set_data(quadratureFieldPoints);
                continue; // Skip the standard per-turn-pair loop for this harmonic
            }
        }

        if (!turnFieldSums) {
            turnFieldSums = build_turn_field_sums();
        }
        const auto& scalePerWinding = turnFieldSums->scalePerHarmonicPerWinding[harmonicIndex];
        const auto& inducingData = inducingFields[harmonicIndex].get_data();

        const auto& inducedDataThisHarmonic = inducedFields[harmonicIndex].get_data();
        for (size_t inducedIndex = 0; inducedIndex < inducedDataThisHarmonic.size(); ++inducedIndex) {
            const auto& inducedFieldPoint = inducedDataThisHarmonic[inducedIndex];
            double totalInducedFieldX = 0;
            double totalInducedFieldY = 0;
            double totalQuadratureInducedFieldX = 0;
            double totalQuadratureInducedFieldY = 0;

            // ROSHEN and SULLIVAN fringing are computed per-point in this loop (not via equivalent current loops)
            // Skip if using ALBACH H-field model since fringing is already added in the ALBACH branch above
            // ALBACH fringing model uses equivalent current loops which are added to inducingFields and processed below,
            // except for gaps beyond its fitted validity, which are routed through the Roshen model here.
            bool albachRoutedGapsPending = _magneticFieldStrengthFringingEffectModel == MagneticFieldStrengthFringingEffectModels::ALBACH &&
                                           !albachOutOfRangeGaps.empty();
            if (!isAlbach && (_magneticFieldStrengthFringingEffectModel == MagneticFieldStrengthFringingEffectModels::ROSHEN ||
                              _magneticFieldStrengthFringingEffectModel == MagneticFieldStrengthFringingEffectModels::SULLIVAN ||
                              albachRoutedGapsPending)) {
                // Every harmonic's gap field, from that harmonic's magnetizing flux (get_gap_source),
                // evaluated only when a functional gap exists to fringe.
                if (includeFringing && hasFunctionalGap && get_gap_source().magnitude > 0) {
                    double magneticFieldStrengthGap = get_gap_source().magnitude;

                    // Multi-column winding: the fringing conventions below are written
                    // for the x>0 window (gaps at x<0 are skipped, edge selection
                    // assumes the point sits right of the gap). Points in the mirrored
                    // (x<0) window see the mirror-symmetric field of the symmetric
                    // gap arrangement: evaluate at the mirrored point and flip the
                    // odd (x) component. Exact for symmetric gapping; asymmetric
                    // lateral gapping across sides is not represented yet.
                    auto fringingInducedPoint = inducedFieldPoint;
                    bool mirroredForFringing = false;
                    if (multiWindowCore && inducedFieldPoint.get_point()[0] < 0) {
                        auto mirroredPoint = inducedFieldPoint.get_point();
                        mirroredPoint[0] = -mirroredPoint[0];
                        fringingInducedPoint.set_point(mirroredPoint);
                        mirroredForFringing = true;
                    }

                    // For ALBACH fringing only the out-of-validity gaps are handled here
                    // (via Roshen); the in-range gaps went through equivalent current loops.
                    auto& gapsToProcess = albachRoutedGapsPending ? albachOutOfRangeGaps : gapping;
                    for (auto& gap : gapsToProcess) {
                        if (gap.get_coordinates().value()[0] < 0) {
                            continue;
                        }
                        // ABT #832: residual (mating-surface) gaps do not contribute
                        // fringing loss -- see the matching gate in the ALBACH branch.
                        if (gap.get_type() != GapType::SUBTRACTIVE && gap.get_type() != GapType::ADDITIVE) {
                            continue;
                        }
                        auto complexFieldPoint = albachRoutedGapsPending ?
                            albachFallbackRoshenModel.get_magnetic_field_strength_between_gap_and_point(gap, magneticFieldStrengthGap, fringingInducedPoint) :
                            _fringingEffectModel->get_magnetic_field_strength_between_gap_and_point(gap, magneticFieldStrengthGap, fringingInducedPoint);

                        // The magnetizing flux's phase: cos in phase, sin in quadrature.
                        double inPhaseFactor = get_gap_source().inPhaseFactor;
                        double quadratureFactor = get_gap_source().quadratureFactor;
                        double fringingX = mirroredForFringing ? -complexFieldPoint.get_real() : complexFieldPoint.get_real();
                        totalInducedFieldX += fringingX * inPhaseFactor;
                        totalInducedFieldY += complexFieldPoint.get_imaginary() * inPhaseFactor;
                        totalQuadratureInducedFieldX += fringingX * quadratureFactor;
                        totalQuadratureInducedFieldY += complexFieldPoint.get_imaginary() * quadratureFactor;
                        if (std::isnan(complexFieldPoint.get_real())) {
                            throw NaNResultException("NaN found in fringing field calculation");
                        }
                        if (std::isnan(complexFieldPoint.get_imaginary())) {
                            throw NaNResultException("NaN found in fringing field calculation");
                        }
                    }
                }
            }

            // Same fringing-only rule as the ALBACH branch above (ABT #182): width
            // samples must not receive turn fields — filament representations of
            // neighbouring flat conductors misread as kA/m of perpendicular field.
            // Fringing contributions (added into totalInducedFieldX/Y above, and the
            // Albach equivalent-current loops, which carry no turn_index) still apply.
            // The turns' field comes from turnFieldSums (above), which applies this rule,
            // the multi-column screening and the own-turn and inside-core exclusions.
            for (size_t windingIndex = 0; windingIndex < scalePerWinding.size(); ++windingIndex) {
                auto [windingFieldX, windingFieldY] = turnFieldSums->fieldPerInducedPointPerWinding[inducedIndex][windingIndex];
                // A turn of winding k carries value * exp(j phase_k): in phase and in quadrature.
                double scale = scalePerWinding[windingIndex];
                totalInducedFieldX += windingFieldX * scale * inPhaseFactorPerWinding[windingIndex];
                totalInducedFieldY += windingFieldY * scale * inPhaseFactorPerWinding[windingIndex];
                totalQuadratureInducedFieldX += windingFieldX * scale * quadratureFactorPerWinding[windingIndex];
                totalQuadratureInducedFieldY += windingFieldY * scale * quadratureFactorPerWinding[windingIndex];
            }

            // Sources that are not turns -- equivalent fringing sources -- are summed per harmonic.
            for (const auto& inducingFieldPoint : inducingData) {
                if (inducingFieldPoint.get_turn_index()) {
                    continue;
                }
                // Multi-column winding: the main column magnetically screens the two
                // window sides from each other (the mirror-image walls).
                if (multiWindowCore &&
                    inducingFieldPoint.get_point()[0] * inducedFieldPoint.get_point()[0] < 0) {
                    continue;
                }

                auto [inducedFieldX, inducedFieldY] = _model->get_magnetic_field_strength_components_between_two_points(inducingFieldPoint, inducedFieldPoint, std::nullopt);

                // An equivalent fringing source: the magnetizing field, on the magnetizing
                // flux's phase.
                double inPhaseFactor = get_gap_source().inPhaseFactor;
                double quadratureFactor = get_gap_source().quadratureFactor;
                totalInducedFieldX += inducedFieldX * inPhaseFactor;
                totalInducedFieldY += inducedFieldY * inPhaseFactor;
                totalQuadratureInducedFieldX += inducedFieldX * quadratureFactor;
                totalQuadratureInducedFieldY += inducedFieldY * quadratureFactor;
                if (std::isnan(inducedFieldX)) {
                    throw NaNResultException("NaN found in magnetic field calculation");
                }
                if (std::isnan(inducedFieldY)) {
                    throw NaNResultException("NaN found in magnetic field calculation");
                }
            }
            ComplexFieldPoint complexFieldPoint;
            complexFieldPoint.set_point(inducedFieldPoint.get_point());
            complexFieldPoint.set_real(totalInducedFieldX);
            complexFieldPoint.set_imaginary(totalInducedFieldY);
            if (inducedFieldPoint.get_turn_index()) {
                complexFieldPoint.set_turn_index(inducedFieldPoint.get_turn_index().value());
            }
            if (inducedFieldPoint.get_label()) {
                complexFieldPoint.set_label(inducedFieldPoint.get_label().value());
            }
            fieldPoints.push_back(complexFieldPoint);

            ComplexFieldPoint quadratureFieldPoint(complexFieldPoint);
            quadratureFieldPoint.set_real(totalQuadratureInducedFieldX);
            quadratureFieldPoint.set_imaginary(totalQuadratureInducedFieldY);
            quadratureFieldPoints.push_back(quadratureFieldPoint);
        }
        complexFieldPerHarmonic[harmonicIndex].set_data(fieldPoints);
        quadratureComplexFieldPerHarmonic[harmonicIndex].set_data(quadratureFieldPoints);
    }

    WindingWindowMagneticStrengthFieldPhasorOutput windingWindowMagneticStrengthFieldOutput;
    windingWindowMagneticStrengthFieldOutput.set_field_per_frequency(complexFieldPerHarmonic);
    windingWindowMagneticStrengthFieldOutput.set_quadrature_field_per_frequency(quadratureComplexFieldPerHarmonic);
    windingWindowMagneticStrengthFieldOutput.set_method_used(to_string(_magneticFieldStrengthModel));
    windingWindowMagneticStrengthFieldOutput.set_origin(ResultOrigin::SIMULATION);
    return windingWindowMagneticStrengthFieldOutput;
}

ComplexFieldPoint MagneticFieldStrengthWangModel::get_magnetic_field_strength_between_two_points(const FieldPoint& inducingFieldPoint, const FieldPoint& inducedFieldPoint, std::optional<size_t> inducingWireIndex) {
    double Hx = 0;
    double Hy = 0;
    if (!inducingWireIndex) {
        return MagneticFieldStrengthLammeranerModel().get_magnetic_field_strength_between_two_points(inducingFieldPoint, inducedFieldPoint);
    }
    else {
        double c = 0;
        double h = 0;
        if (_wirePerWinding[inducingWireIndex.value()].get_type() == WireType::FOIL) {
            c = resolve_dimensional_values(_wirePerWinding[inducingWireIndex.value()].get_conducting_width().value());
            h = resolve_dimensional_values(_wirePerWinding[inducingWireIndex.value()].get_conducting_height().value());
        }
        else {
            h = resolve_dimensional_values(_wirePerWinding[inducingWireIndex.value()].get_conducting_width().value());
            c = resolve_dimensional_values(_wirePerWinding[inducingWireIndex.value()].get_conducting_height().value());
        }
        double k = c / h;
        double lambda = 0.01 * k + 0.66;

        if (inducedFieldPoint.get_label() && inducingFieldPoint.get_label()) {
            auto inducingLabel = inducingFieldPoint.get_label().value();
            auto inducedLabel = inducedFieldPoint.get_label().value();
            auto current = inducingFieldPoint.get_value();
            double distanceX = inducingFieldPoint.get_point()[0] - inducedFieldPoint.get_point()[0];
            double distanceY = inducingFieldPoint.get_point()[1] - inducedFieldPoint.get_point()[1];
            double distance = hypot(distanceX, distanceY);
            if (inducedLabel == "widthsample") {
                // Width samples are generic spatial points, not Wang's paired edge
                // points, so the label-pair formulas below do not apply. Treat each
                // Wang inducing point as a filament via Lammeraner. Wang's inducing
                // points carry the full turn current (its own formulas halve it), so
                // halve the filament current here to keep the total at I per turn.
                FieldPoint halfCurrentPoint = inducingFieldPoint;
                halfCurrentPoint.set_value(current * 0.5);
                return MagneticFieldStrengthLammeranerModel().get_magnetic_field_strength_between_two_points(halfCurrentPoint, inducedFieldPoint);
            }
            if (inducingLabel == "left") {
                if (inducedLabel == "left") {
                    double tetha1 = asin(fabs(distanceY) / distance);
                    Hy = 0.5 * current / (2 * std::numbers::pi * lambda * h) + 0.5 * current * cos(tetha1) / (2 * std::numbers::pi * sqrt(pow(lambda * h, 2) + pow(distanceY, 2)));
                }
                else if (inducedLabel == "right") {
                    double tetha2 = asin(fabs(distanceY) / distance);
                    Hy = -0.5 * current / (2 * std::numbers::pi * (c - lambda * h)) - 0.5 * current * cos(tetha2) / (2 * std::numbers::pi * sqrt(pow(c - lambda * h, 2) + pow(distanceY, 2)));
                }
                else if (inducedLabel == "top") {
                    if (distanceY > 0) {
                        Hx = 0;
                    }
                    else {
                        Hx = (current - Hy * 2 * h) / (2 * c);
                    }
                }
                else if (inducedLabel == "bottom") {
                    if (distanceY > 0) {
                        Hx = (current - Hy * 2 * h) / (2 * c);
                    }
                    else {
                        Hx = 0;
                    }
                }
                else {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT, "Wrong inducedLabel: " + inducedLabel);
                }
            }
            else if (inducingLabel == "right") {
                if (inducedLabel == "right") {
                    double tetha1 = asin(fabs(distanceY) / distance);
                    Hy = 0.5 * current / (2 * std::numbers::pi * lambda * h) + 0.5 * current * cos(tetha1) / (2 * std::numbers::pi * sqrt(pow(lambda * h, 2) + pow(distanceY, 2)));
                }
                else if (inducedLabel == "left") {
                    double tetha2 = asin(fabs(distanceY) / distance);
                    Hy = -0.5 * current / (2 * std::numbers::pi * (c - lambda * h)) - 0.5 * current * cos(tetha2) / (2 * std::numbers::pi * sqrt(pow(c - lambda * h, 2) + pow(distanceY, 2)));
                }
                else if (inducedLabel == "bottom") {
                    if (distanceY > 0) {
                        Hx = 0;
                    }
                    else {
                        Hx = (current - Hy * 2 * h) / (2 * c);
                    }
                }
                else if (inducedLabel == "top") {
                    if (distanceY > 0) {
                        Hx = (current - Hy * 2 * h) / (2 * c);
                    }
                    else {
                        Hx = 0;
                    }
                }
            }
            else if (inducingLabel == "bottom") {
                if (inducedLabel == "bottom") {
                    double tetha1 = asin(fabs(distanceY) / distance);
                    Hy = 0.5 * current / (2 * std::numbers::pi * lambda * h) + 0.5 * current * cos(tetha1) / (2 * std::numbers::pi * sqrt(pow(lambda * h, 2) + pow(distanceY, 2)));
                }
                else if (inducedLabel == "top") {
                    double tetha2 = asin(fabs(distanceY) / distance);
                    Hy = -0.5 * current / (2 * std::numbers::pi * (c - lambda * h)) - 0.5 * current * cos(tetha2) / (2 * std::numbers::pi * sqrt(pow(c - lambda * h, 2) + pow(distanceY, 2)));
                }
                else if (inducedLabel == "right") {
                    if (distanceY > 0) {
                        Hx = 0;
                    }
                    else {
                        Hx = (current - Hy * 2 * h) / (2 * c);
                    }
                }
                else if (inducedLabel == "left") {
                    if (distanceY > 0) {
                        Hx = (current - Hy * 2 * h) / (2 * c);
                    }
                    else {
                        Hx = 0;
                    }
                }
                else {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT, "Wrong inducedLabel: " + inducedLabel);
                }
            }
            else if (inducingLabel == "top") {
                if (inducedLabel == "top") {
                    double tetha1 = asin(fabs(distanceY) / distance);
                    Hy = 0.5 * current / (2 * std::numbers::pi * lambda * h) + 0.5 * current * cos(tetha1) / (2 * std::numbers::pi * sqrt(pow(lambda * h, 2) + pow(distanceY, 2)));
                }
                else if (inducedLabel == "bottom") {
                    double tetha2 = asin(fabs(distanceY) / distance);
                    Hy = -0.5 * current / (2 * std::numbers::pi * (c - lambda * h)) - 0.5 * current * cos(tetha2) / (2 * std::numbers::pi * sqrt(pow(c - lambda * h, 2) + pow(distanceY, 2)));
                }
                else if (inducedLabel == "left") {
                    if (distanceY > 0) {
                        Hx = 0;
                    }
                    else {
                        Hx = (current - Hy * 2 * h) / (2 * c);
                    }
                }
                else if (inducedLabel == "right") {
                    if (distanceY > 0) {
                        Hx = (current - Hy * 2 * h) / (2 * c);
                    }
                    else {
                        Hx = 0;
                    }
                }
            }
            else {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, "Wrong inducingLabel: " + inducingLabel);
            }
        }
        else {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Wang Magnetic Field model must be used with his CoilMesher model");
        }
    }

    ComplexFieldPoint complexFieldPoint;
    complexFieldPoint.set_imaginary(Hy);
    complexFieldPoint.set_point(inducedFieldPoint.get_point());
    complexFieldPoint.set_real(Hx);
    if (inducedFieldPoint.get_turn_index()) {
        complexFieldPoint.set_turn_index(inducedFieldPoint.get_turn_index().value());
    }
    if (inducedFieldPoint.get_turn_length()) {
        complexFieldPoint.set_turn_length(inducedFieldPoint.get_turn_length().value());
    }
    return complexFieldPoint;   
}

std::pair<double, double> MagneticFieldStrengthBinnsLawrensonModel::get_magnetic_field_strength_components_between_two_points(const FieldPoint& inducingFieldPoint, const FieldPoint& inducedFieldPoint, std::optional<size_t> inducingWireIndex) {
    double Hx;
    double Hy;

    double distanceX = inducingFieldPoint.get_point()[0] - inducedFieldPoint.get_point()[0];
    double distanceY = inducingFieldPoint.get_point()[1] - inducedFieldPoint.get_point()[1];

    bool computeRectangularFiled = false;

    if (inducingWireIndex) {
        if (_wirePerWinding[inducingWireIndex.value()].get_type() != WireType::ROUND && _wirePerWinding[inducingWireIndex.value()].get_type() != WireType::LITZ) {
            computeRectangularFiled = true;
        }
        else {
            double wireRadius = _wireMaxOuterWidth[inducingWireIndex.value()] / 2;
            // Inside the wire the current enclosed at radius r is I r^2 / R^2 (uniform density), so the
            // field is I r / (2 pi R^2). It was zeroed, which dropped the conductor's own stored energy from
            // the leakage integral (ABT #1240).
            double divisor = 2 * std::numbers::pi * std::max(pow(distanceY, 2) + pow(distanceX, 2), wireRadius * wireRadius);
            Hx = -inducingFieldPoint.get_value() * (distanceY) / divisor;
            Hy = inducingFieldPoint.get_value() * (distanceX) / divisor;
            if (std::isnan(Hx) || std::isnan(Hy)) {
                throw NaNResultException("NaN found in Binns Lawrenson's model for magnetic field");
            }
        }
    }
    else {
        double divisor = 2 * std::numbers::pi * (pow(distanceY, 2) + pow(distanceX, 2));
        Hx = -inducingFieldPoint.get_value() * (distanceY) / divisor;
        Hy = inducingFieldPoint.get_value() * (distanceX) / divisor;
        if (std::isnan(Hx) || std::isnan(Hy)) {
            throw NaNResultException("NaN found in Binns Lawrenson's model for magnetic field");
        }
    }


    if (computeRectangularFiled) {
        double a = resolve_dimensional_values(_wirePerWinding[inducingWireIndex.value()].get_conducting_width().value()) / 2;
        double b = resolve_dimensional_values(_wirePerWinding[inducingWireIndex.value()].get_conducting_height().value()) / 2;
        double x = inducedFieldPoint.get_point()[0] - inducingFieldPoint.get_point()[0];
        double y = inducedFieldPoint.get_point()[1] - inducingFieldPoint.get_point()[1];

        if (inducingFieldPoint.get_rotation()) {
            double modulo = hypot(x, y);
            double currentAngle = atan2(y, x);
            double turnAngle = inducingFieldPoint.get_rotation().value() / 180 * std::numbers::pi;
            if (currentAngle < 0) {
                currentAngle += 2 * std::numbers::pi;
            }
            double totalAngle = currentAngle - turnAngle;

            x = modulo * cos(totalAngle);
            y = modulo * sin(totalAngle);
        }

        // Field of a uniform current I over the rectangle |x'| < a, |y'| < b, from the double integral of the
        // line-current kernel. With u = x - x', v = y - y' and P(u, v) = u ln(u^2 + v^2) / 2 + v atan(u / v), a
        // primitive of v / (u^2 + v^2) du dv, Hx = J/(2 pi) S(u, v) and Hy = -J/(2 pi) S(v, u) in this model's sign
        // convention (the filament branch above), S being P summed over the four corners. It is exact inside the
        // conductor too. The previous atan/log form with quadrant corrections (identical outside, to 1e-6 on 2e4
        // random points) zeroed the field in the copper, which dropped the tracks' own energy from the planar
        // leakage integral (ABT #1240).
        auto primitive = [](double u, double v) {
            double squaredRadius = u * u + v * v;
            double logarithmTerm = squaredRadius > 0 ? 0.5 * u * log(squaredRadius) : 0.0;
            double arctangentTerm = v != 0 ? v * atan(u / v) : 0.0;
            return logarithmTerm + arctangentTerm;
        };
        auto cornerSum = [&primitive](double u1, double u2, double v1, double v2) {
            return primitive(u2, v2) - primitive(u1, v2) - primitive(u2, v1) + primitive(u1, v1);
        };
        double currentDensityOver2Pi = inducingFieldPoint.get_value() / (4.0 * a * b) / (2.0 * std::numbers::pi);
        Hx = currentDensityOver2Pi * cornerSum(x - a, x + a, y - b, y + b);
        Hy = -currentDensityOver2Pi * cornerSum(y - b, y + b, x - a, x + a);
        if (std::isnan(Hx) || std::isnan(Hy)) {
            throw NaNResultException("NaN found in Binns Lawrenson's model for magnetic field");
        }
    }


    if (inducingFieldPoint.get_rotation()) {
        double modulo = hypot(Hx, Hy);
        double currentAngle = atan2(Hy, Hx);
        if (currentAngle < 0) {
            currentAngle += 2 * std::numbers::pi;
        }
        double turnAngle = inducingFieldPoint.get_rotation().value() / 180 * std::numbers::pi;
        double totalAngle = currentAngle + turnAngle;
        Hx = modulo * cos(totalAngle);
        Hy = modulo * sin(totalAngle);
        if (std::isnan(Hx) || std::isnan(Hy)) {
            throw NaNResultException("NaN found in Binns Lawrenson's model for magnetic field");
        }
    }

    return {Hx, Hy};
}

ComplexFieldPoint MagneticFieldStrengthBinnsLawrensonModel::get_magnetic_field_strength_between_two_points(const FieldPoint& inducingFieldPoint, const FieldPoint& inducedFieldPoint, std::optional<size_t> inducingWireIndex) {
    auto [Hx, Hy] = get_magnetic_field_strength_components_between_two_points(inducingFieldPoint, inducedFieldPoint, inducingWireIndex);
    ComplexFieldPoint complexFieldPoint;
    complexFieldPoint.set_imaginary(Hy);
    complexFieldPoint.set_point(inducedFieldPoint.get_point());
    complexFieldPoint.set_real(Hx);
    if (inducedFieldPoint.get_turn_index()) {
        complexFieldPoint.set_turn_index(inducedFieldPoint.get_turn_index().value());
    }
    if (inducedFieldPoint.get_turn_length()) {
        complexFieldPoint.set_turn_length(inducedFieldPoint.get_turn_length().value());
    }
    return complexFieldPoint;   
}

// ABT #376: Dowell's one-dimensional field (see the class comment in the header). The layers are
// assumed to span the winding breadth b, so the field is PARALLEL to them and depends only on the
// ampere-turns enclosed between the induced point and the zero-field boundary. Per inducing
// conductor that is a step: the conductor's own I/b is felt on one side of it and nothing on the
// other, so summing over the conductors of a layered winding rebuilds Dowell's MMF staircase.
//
// Orientation: MKF's concentric windows stack layers along x (radial) with the field running
// along y (axial), which is Dowell's own arrangement, so the step is taken on the x coordinate
// and the field is returned on y. Toroidal/other layouts are not Dowell's geometry and the
// caller should choose a two-dimensional model there.
ComplexFieldPoint MagneticFieldStrengthDowellModel::get_magnetic_field_strength_between_two_points(const FieldPoint& inducingFieldPoint, const FieldPoint& inducedFieldPoint, std::optional<size_t> inducingWireIndex) {
    ComplexFieldPoint magneticFieldStrengthPoint;
    magneticFieldStrengthPoint.set_point(inducedFieldPoint.get_point());
    if (inducedFieldPoint.get_label()) {
        magneticFieldStrengthPoint.set_label(inducedFieldPoint.get_label().value());
    }

    // Without a breadth there is no Dowell field to speak of; refusing beats inventing one.
    if (_windingWindowBreadth <= 0) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "Dowell's field model needs the winding breadth (the window dimension parallel to the "
            "layers); none was supplied, so H = MMF / b is undefined");
    }

    double inducingRadialPosition = inducingFieldPoint.get_point()[0];
    double inducedRadialPosition = inducedFieldPoint.get_point()[0];
    double fieldStep = inducingFieldPoint.get_value() / _windingWindowBreadth;

    // The conductor contributes to points OUTSIDE it (further from the zero-field boundary) and
    // not to points inside. On the conductor itself the enclosed fraction is taken as half, which
    // is Dowell's own treatment of the layer carrying the current.
    double contribution;
    if (inducedRadialPosition > inducingRadialPosition) {
        contribution = fieldStep;
    }
    else if (inducedRadialPosition < inducingRadialPosition) {
        contribution = 0;
    }
    else {
        contribution = fieldStep / 2;
    }

    magneticFieldStrengthPoint.set_real(0);
    magneticFieldStrengthPoint.set_imaginary(contribution);
    return magneticFieldStrengthPoint;
}

std::pair<double, double> MagneticFieldStrengthLammeranerModel::get_magnetic_field_strength_components_between_two_points(const FieldPoint& inducingFieldPoint, const FieldPoint& inducedFieldPoint, std::optional<size_t> inducingWireIndex) {
    if (inducingWireIndex && _wirePerWinding[inducingWireIndex.value()].get_type() != WireType::ROUND && _wirePerWinding[inducingWireIndex.value()].get_type() != WireType::LITZ) {
        auto complexFieldPoint = MagneticFieldStrengthBinnsLawrensonModel().get_magnetic_field_strength_between_two_points(inducingFieldPoint, inducedFieldPoint);
        return {complexFieldPoint.get_real(), complexFieldPoint.get_imaginary()};
    }

    double turnLength = 1;
    if (inducingFieldPoint.get_turn_length()) {
        turnLength = inducingFieldPoint.get_turn_length().value();
    }
    const auto& inducedPoint = inducedFieldPoint.get_point();
    const auto& inducingPoint = inducingFieldPoint.get_point();
    double distanceX = inducedPoint[0] - inducingPoint[0];
    double distanceY = inducedPoint[1] - inducingPoint[1];
    double distance = hypot(distanceY, distanceX);

    if (inducingWireIndex && distance < _wireMaxOuterWidth[inducingWireIndex.value()] / 2) {
        return {0, 0};
    }

    // The azimuthal field of a line current. Its unit vector is (cos(angle + pi/2),
    // sin(angle + pi/2)) with angle = atan2(dy, dx), which is (-dy, dx) / distance: the same
    // direction without an atan2, a cos and a sin per point pair. The leading minus in the
    // module gives the same clockwise field a line current has in the Binns-Lawrenson model.
    double magneticFiledStrengthModule = -inducingFieldPoint.get_value() / 2 / std::numbers::pi / distance * turnLength / hypot(turnLength, distance);
    double Hx;
    double Hy;
    if (distance > 0) {
        Hx = magneticFiledStrengthModule * (-distanceY / distance);
        Hy = magneticFiledStrengthModule * (distanceX / distance);
    }
    else {
        // Coincident points (only reachable for a non-turn source): angle = atan2(0, 0) = 0.
        Hx = magneticFiledStrengthModule * cos(std::numbers::pi / 2);
        Hy = magneticFiledStrengthModule * sin(std::numbers::pi / 2);
    }

    if (std::isnan(Hx) || std::isnan(Hy)) {
        throw NaNResultException("NaN found in Lammeraner's model for magnetic field");
    }
    return {Hx, Hy};
}

ComplexFieldPoint MagneticFieldStrengthLammeranerModel::get_magnetic_field_strength_between_two_points(const FieldPoint& inducingFieldPoint, const FieldPoint& inducedFieldPoint, std::optional<size_t> inducingWireIndex) {
    auto [Hx, Hy] = get_magnetic_field_strength_components_between_two_points(inducingFieldPoint, inducedFieldPoint, inducingWireIndex);

    ComplexFieldPoint complexFieldPoint;
    complexFieldPoint.set_imaginary(Hy);
    complexFieldPoint.set_point(inducedFieldPoint.get_point());
    complexFieldPoint.set_real(Hx);
    if (inducedFieldPoint.get_turn_index()) {
        complexFieldPoint.set_turn_index(inducedFieldPoint.get_turn_index().value());
    }
    if (inducedFieldPoint.get_turn_length()) {
        complexFieldPoint.set_turn_length(inducedFieldPoint.get_turn_length().value());
    }
    return complexFieldPoint;
}

bool MagneticFieldStrengthAlbachModel::is_gap_within_validity_range(CoreGap gap) {
    if (!gap.get_section_dimensions() || !gap.get_coordinates()) {
        return false;
    }
    double rc = gap.get_section_dimensions().value()[0] / 2;
    double xi = gap.get_length() / (2 * rc);
    double x = 1 - 1.05 * xi - 2.88 * pow(xi, 2) - 8.8 * pow(xi, 3);
    double denominator = 0.25 - 1.569 * xi + 4.34 * pow(xi, 2) - 7.042 * pow(xi, 3);
    return x >= 0 && denominator > 0;
}

FieldPoint MagneticFieldStrengthAlbachModel::get_equivalent_inducing_point_for_gap(CoreGap gap, double magneticFieldStrengthGap) {
    if (!gap.get_section_dimensions()) {
        throw GapException("Gap is missing section dimensions");
    }
    if (!gap.get_coordinates()) {
        throw GapException("Gap is missing coordinates");
    }
    double rc = gap.get_section_dimensions().value()[0] / 2;
    double xi = gap.get_length() / (2 * rc);
    double x = 1 - 1.05 * xi - 2.88 * pow(xi, 2) - 8.8 * pow(xi, 3);
    if (x < 0) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something went wrong with Albach method with x");
    }
    // The fitted denominator polynomial crosses zero at xi ~ 0.2755 — BEFORE the x < 0
    // guard above fires (xi ~ 0.42) — so gaps in that band produced a diverging,
    // sign-flipping equivalent fringing current with no error. Albach's fit is only
    // valid for small xi; fail loudly outside it instead of returning garbage.
    double denominator = 0.25 - 1.569 * xi + 4.34 * pow(xi, 2) - 7.042 * pow(xi, 3);
    if (denominator <= 0) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Albach fringing model out of its validity range: gap length / column diameter ratio too large (xi = " + std::to_string(xi) + ")");
    }
    double current = (magneticFieldStrengthGap * gap.get_length()) / denominator;
    double eta = x * rc;

    if (eta > gap.get_section_dimensions().value()[0] / 2) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something went wrong with Albach method with eta");
    }
    FieldPoint fieldPoint;
    // Position the equivalent wire at distance eta from the center leg axis
    // For center leg gaps (x=0), the wire is placed at radius eta
    // For lateral gaps (x>0 or x<0), adjust based on gap position
    if (gap.get_coordinates().value()[0] > 0) {
        // Gap on positive x side - wire is at gap_x - eta (closer to center)
        fieldPoint.set_point({gap.get_coordinates().value()[0] - eta, gap.get_coordinates().value()[1]});
    }
    else if (gap.get_coordinates().value()[0] < 0) {
        // Gap on negative x side - wire is at gap_x + eta (closer to center)
        fieldPoint.set_point({gap.get_coordinates().value()[0] + eta, gap.get_coordinates().value()[1]});
    }
    else {
        // Center leg gap (x=0) - wire is placed at radius eta from axis
        fieldPoint.set_point({eta, gap.get_coordinates().value()[1]});
    }
    fieldPoint.set_value(current);
    return fieldPoint;
}

// SIGN (ABT #1463). The field returned is the gap's fringing field in the frame of the turn
// kernels, for a POSITIVE magnetizing MMF (into-dot current): flux up (+y) the central leg and
// down (-y) the lateral legs. Every turn kernel (Albach air coil, Binns-Lawrenson, Lammeraner)
// gives +y at the central leg for a positive filament current, so a positive MMF must give a gap
// whose flux runs +y there, and its fringing field is the continuation of that flux.
// Roshen's conformal expression below, taken as written with H_g > 0 and xi the distance from the
// leg surface into the window, reads at the gap mouth (xi -> 0+, dy = 0): Hy = -0.9 H_g, i.e.
// ANTIPARALLEL to the gap's own field; and above the gap centre (dy > 0) its Hx points away from
// the leg, where the flux of a +y gap bends back into the upper pole face (Hx towards the leg). So
// as written it is the field of a gap whose flux runs -y: it is negated for a +y gap. With the old
// sign the fringing field subtracted from the turns' field where the two add, and the cross term
// had the wrong sign: on a synthetic DCM flyback (ETD29, 0.82 mm gap, magnetizing MMF = the
// windings' true net MMF) h = 1 R_ac/R_dc read primary 1.48x / secondary 3.16x OMFEM with the old
// sign and 0.88x / 0.95x with this one.
// A lateral leg's window lies on its -x side and its flux runs -y: the central-leg field is
// mirrored (x -> -x flips Hx) and reversed (both components flip), so Hx keeps the central sign
// and Hy flips. The old code fed the lateral dx < 0 straight into the central expression, whose
// atan branch (m) is built for dx > 0: its Hy jumped from -1.35 H_g to +0.45 H_g across the circle
// of radius l_g/2 around the edge.
ComplexFieldPoint MagneticFieldStrengthRoshenModel::get_magnetic_field_strength_between_gap_and_point(CoreGap gap, double magneticFieldStrengthGap, FieldPoint inducedFieldPoint) {
    bool centralGap = gap.get_coordinates().value()[0] == 0;
    double distanceFromCenterEdgeGapX;
    if (centralGap) {
        distanceFromCenterEdgeGapX = inducedFieldPoint.get_point()[0] - (gap.get_coordinates().value()[0] + gap.get_section_dimensions().value()[0] / 2);
    }
    else {
        // Measured from the lateral leg's window-side surface INTO the window (towards -x).
        distanceFromCenterEdgeGapX = (gap.get_coordinates().value()[0] - gap.get_section_dimensions().value()[0] / 2) - inducedFieldPoint.get_point()[0];
    }
    double distanceFromCenterEdgeGapY = inducedFieldPoint.get_point()[1] - gap.get_coordinates().value()[1];
    double halfGapLength = gap.get_length() / 2;

    double magneticIntensityXDividend = pow(distanceFromCenterEdgeGapX, 2) + pow(distanceFromCenterEdgeGapY - halfGapLength, 2);
    double magneticIntensityXDivisor = pow(distanceFromCenterEdgeGapX, 2) + pow(distanceFromCenterEdgeGapY + halfGapLength, 2);
    double Hx = -0.9 * magneticFieldStrengthGap / 2 / std::numbers::pi * log(magneticIntensityXDividend / magneticIntensityXDivisor);

    double m;
    if (pow(distanceFromCenterEdgeGapX, 2) + pow(distanceFromCenterEdgeGapY, 2) > pow(halfGapLength, 2)) {
        m = 0;
    }
    else {
        m = 1;
    }

    double x = distanceFromCenterEdgeGapX * halfGapLength / (pow(distanceFromCenterEdgeGapX, 2) + pow(distanceFromCenterEdgeGapY, 2) - pow(halfGapLength, 2));
    double Hy = -0.9 * magneticFieldStrengthGap / std::numbers::pi * (atan(x) + m * std::numbers::pi);

    // (Hx, Hy) as written: a gap with flux along -y and its window at +xi (see SIGN above).
    // Central leg, flux +y, window at +x: the negative, (-Hx, -Hy).
    // Lateral leg, flux -y, window at -x: the central field mirrored (-Hx, -Hy) -> (Hx, -Hy),
    // then reversed (flux -y) -> (-Hx, Hy).
    Hx = -Hx;
    if (centralGap) {
        Hy = -Hy;
    }

    ComplexFieldPoint complexFieldPoint;
    complexFieldPoint.set_imaginary(Hy);
    complexFieldPoint.set_point(inducedFieldPoint.get_point());
    complexFieldPoint.set_real(Hx);
    if (inducedFieldPoint.get_turn_index()) {
        complexFieldPoint.set_turn_index(inducedFieldPoint.get_turn_index().value());
    }
    return complexFieldPoint;
}


// ============================================================================
// MagneticFieldStrengthAlbach2DModel Implementation (Air Coil / Biot-Savart)
// ============================================================================

ComplexFieldPoint MagneticFieldStrengthAlbach2DModel::get_magnetic_field_strength_between_two_points(
    const FieldPoint& inducingFieldPoint,
    const FieldPoint& inducedFieldPoint, 
    std::optional<size_t> inducingWireIndex
) {
    // ALBACH model calculates field from all turns at once via calculateTotalFieldAtPoint()
    // This per-turn-pair method should not be called
    (void)inducingFieldPoint;
    (void)inducedFieldPoint;
    (void)inducingWireIndex;
    throw std::runtime_error("ALBACH model does not support per-turn-pair field calculation. Use calculateTotalFieldAtPoint() instead.");
}

std::pair<double, double> MagneticFieldStrengthAlbach2DModel::calculateMagneticField(double r, double z) const {
    // Calculate H directly using analytical Biot-Savart formulas for circular filaments
    // Uses complete elliptic integrals of the first and second kind
    
    // Handle r near zero to avoid division by zero
    if (r < 1e-10) {
        return {0.0, 0.0};
    }
    
    double H_r_total = 0.0;
    double H_z_total = 0.0;
    
    for (const auto& turn : _turns) {
        double I = turn.current;
        if (std::abs(I) < 1e-15) continue;
        
        if (turn.isRectangular()) {
            // Rectangular conductor: use filamentary subdivision
            double width = turn.width;
            double height = turn.height;
            // PERF-001: Adaptive subdivision (3-8) based on dim/skinDepth
            int numR = std::max(3, std::min(8, static_cast<int>(std::ceil(width / turn.skinDepth))));
            int numZ = std::max(3, std::min(8, static_cast<int>(std::ceil(height / turn.skinDepth))));
            double dI = I / (numR * numZ);
            
            for (int ir = 0; ir < numR; ++ir) {
                for (int iz = 0; iz < numZ; ++iz) {
                    double fr = (ir + 0.5) / numR;
                    double fz = (iz + 0.5) / numZ;
                    
                    double rf = turn.r - width/2 + width * fr;
                    double zf = turn.z - height/2 + height * fz;
                    
                    if (rf < 1e-10) continue;
                    
                    // Direct Biot-Savart calculation for circular filament
                    double deltaZ = z - zf;
                    double sumR = r + rf;
                    double diffR = r - rf;
                    
                    double denom = sumR * sumR + deltaZ * deltaZ;
                    if (denom < 1e-20) continue;
                    
                    double k2 = 4 * r * rf / denom;
                    double k = std::sqrt(k2);
                    if (k > 0.999999) k = 0.999999;
                    
                    if (k > 1e-10) {
                        double K_k = std::comp_ellint_1(k);
                        double E_k = std::comp_ellint_2(k);
                        
                        double sqrtDenom = std::sqrt(denom);
                        double denomDiffR = diffR * diffR + deltaZ * deltaZ;
                        
                        if (denomDiffR > 1e-20) {
                            double prefactor = dI / (2 * std::numbers::pi);
                            
                            H_r_total += prefactor * deltaZ / (r * sqrtDenom) * 
                                  (-K_k + E_k * (rf*rf + r*r + deltaZ*deltaZ) / denomDiffR);
                            
                            H_z_total += prefactor / sqrtDenom * 
                                  (K_k + E_k * (rf*rf - r*r - deltaZ*deltaZ) / denomDiffR);
                        }
                    }
                }
            }
        } else {
            // Round wire: single filament
            double r0 = turn.r;
            double z0 = turn.z;
            
            double deltaZ = z - z0;
            double sumR = r + r0;
            double diffR = r - r0;
            
            double denom = sumR * sumR + deltaZ * deltaZ;
            if (denom < 1e-20) continue;
            
            double k2 = 4 * r * r0 / denom;
            double k = std::sqrt(k2);
            if (k > 0.999999) k = 0.999999;
            
            if (k > 1e-10) {
                double K_k = std::comp_ellint_1(k);
                double E_k = std::comp_ellint_2(k);
                
                double sqrtDenom = std::sqrt(denom);
                double denomDiffR = diffR * diffR + deltaZ * deltaZ;
                
                if (denomDiffR > 1e-20) {
                    double prefactor = I / (2 * std::numbers::pi);
                    
                    H_r_total += prefactor * deltaZ / (r * sqrtDenom) * 
                          (-K_k + E_k * (r0*r0 + r*r + deltaZ*deltaZ) / denomDiffR);
                    
                    H_z_total += prefactor / sqrtDenom * 
                          (K_k + E_k * (r0*r0 - r*r - deltaZ*deltaZ) / denomDiffR);
                }
            }
        }
    }
    
    return {H_r_total, H_z_total};
}

ComplexFieldPoint MagneticFieldStrengthAlbach2DModel::calculateTotalFieldAtPoint(FieldPoint inducedFieldPoint) {
    // Extract induced point coordinates - in 2D cross section: [0] = x (radial), [1] = y (axial)
    double r = std::abs(inducedFieldPoint.get_point()[0]);
    double z = inducedFieldPoint.get_point()[1];
    
    // Calculate H field from all turns
    auto [H_r, H_z] = calculateMagneticField(r, z);
    
    // Convert to 2D Cartesian: real = radial (Hx), imaginary = axial (Hy)
    ComplexFieldPoint result;
    result.set_real(H_r);
    result.set_imaginary(H_z);
    result.set_point(inducedFieldPoint.get_point());
    if (inducedFieldPoint.get_turn_index()) {
        result.set_turn_index(inducedFieldPoint.get_turn_index().value());
    }
    if (inducedFieldPoint.get_label()) {
        result.set_label(inducedFieldPoint.get_label().value());
    }
    
    return result;
}

void MagneticFieldStrengthAlbach2DModel::setupFromMagnetic(
    Magnetic magnetic, 
    const std::vector<Wire>& wirePerWinding
) {
    // Get turns description
    if (!magnetic.get_coil().get_turns_description()) {
        throw std::runtime_error("Missing turns description in coil");
    }
    auto turns = magnetic.get_coil().get_turns_description().value();
    
    // Set up all turns from the coil
    _turns.clear();
    for (size_t turnIdx = 0; turnIdx < turns.size(); ++turnIdx) {
        auto& turn = turns[turnIdx];
        AlbachTurnPosition albachTurn;
        
        // In 2D cross-section, x = radial, y = axial
        albachTurn.r = std::abs(turn.get_coordinates()[0]);
        albachTurn.z = turn.get_coordinates()[1];
        albachTurn.current = 1.0; // Will be scaled per harmonic
        albachTurn.turnIndex = turnIdx;

        // Get wire info for this turn to set dimensions for rectangular wires
        auto windingIndex = magnetic.get_mutable_coil().get_winding_index_by_name(turn.get_winding());
        if (windingIndex < wirePerWinding.size()) {
            auto& wire = wirePerWinding[windingIndex];
            if (wire.get_type() != WireType::ROUND && wire.get_type() != WireType::LITZ) {
                // Rectangular, foil, or planar wire - set dimensions for subdivision
                if (wire.get_conducting_width()) {
                    albachTurn.width = resolve_dimensional_values(wire.get_conducting_width().value());
                }
                if (wire.get_conducting_height()) {
                    albachTurn.height = resolve_dimensional_values(wire.get_conducting_height().value());
                }
            }
            // For round/litz wires, width and height stay at 0 (point filament)
        }
        _turns.push_back(albachTurn);
    }
}



// ============================================================================
// MagneticFieldStrengthSullivanModel Implementation
// (2D Image Method / Biot-Savart for gap fringing field)
// ============================================================================
//
// THEORY (from Sullivan's shapeopt MATLAB code):
// -------
// The air gap is modeled as a set of current filaments distributed along the
// gap length. For each filament at position R_gap:
//   - A "cross" current (+I_per_div, into page) is placed at the gap face
//   - A "dot" current (-I_per_div, out of page) is placed at the mirror
//     position (reflected about x=0 for center gaps)
//
// The winding window (width bw, height hw) is the fundamental unit cell.
// Image copies are tiled in both x and y:
//   x: at x_center + n * 2*hw,  n in [-imageUnitsX, +imageUnitsX]
//   y: at y_center + m * bw,    m in [-imageUnitsY, +imageUnitsY]
//
// Total field at point P is superposition of all image contributions:
//   B(P) = sum (mu_0*I)/(2*pi) * (P - R_fil) / |P - R_fil|^2
// Then H = B / mu_0
//
// A 0.9 attenuation factor is applied (same as Roshen) to account for
// the fraction of MMF that produces external fringing field.
//
// MAPPING FROM MATLAB CODE:
// -------------------------
// In the original shapeopt code (function Bfinite):
//   - pvec(1:2) = unit_of_X, unit_of_Y  -> _imageUnitsX, _imageUnitsY
//   - pvec(3:4) = bw, hw                -> estimated from gap geometry
//   - pvec(7:8) = gw, gap_div           -> gap.get_length(), _gapDivisions
//   - pvec(9)   = I_per_gap_div         -> I_total / _gapDivisions
//   - pvec(10)  = Rgbase(k)             -> filament position (complex)
//   - const1 = u0*(j)*I/(2*pi)          -> Biot-Savart coefficient
//   - R1_1/R1_2: cross/dot current pair
//   - center_matrix: image unit centers
//

ComplexFieldPoint MagneticFieldStrengthSullivanModel::get_magnetic_field_strength_between_gap_and_point(
    CoreGap gap, double magneticFieldStrengthGap, FieldPoint inducedFieldPoint) {

    if (!gap.get_section_dimensions()) {
        throw GapException("Gap is missing section dimensions");
    }
    if (!gap.get_coordinates()) {
        throw GapException("Gap is missing coordinates");
    }

    // ---- Extract gap geometry ----
    double gapLength = gap.get_length();
    double gapX = gap.get_coordinates().value()[0];
    double gapY = gap.get_coordinates().value()[1];
    double columnWidth = gap.get_section_dimensions().value()[0];

    // ---- Estimate winding window dimensions ----
    // bw: window breadth in y-direction (along the gap), approx = column width
    // hw: window height in x-direction (perpendicular to gap)
    // For center gap (gapX=0): hw ~ columnWidth (the window extends from
    //   the centerpost edge outward)
    // For lateral gap: hw ~ 2*|gapX|
    double bw = columnWidth;
    double hw;
    if (std::abs(gapX) > 1e-10) {
        hw = 2.0 * std::abs(gapX);
    } else {
        hw = columnWidth;
    }

    // ---- Compute total gap current (Ampere's law: NI = H_gap * gapLength) ----
    double I_total = magneticFieldStrengthGap * gapLength;
    double I_per_div = I_total / _gapDivisions;

    // ---- Empirical attenuation (consistent with Roshen) ----
    double attenuationFactor = 0.9;

    // ---- Gap filament spacing ----
    double gapGrid = gapLength / _gapDivisions;

    // ---- Point of interest ----
    double xP = inducedFieldPoint.get_point()[0];
    double yP = inducedFieldPoint.get_point()[1];

    // Accumulate B field components
    double Bx_total = 0.0;
    double By_total = 0.0;

    double u0 = Constants().vacuumPermeability;

    // ---- For each gap filament ----
    for (int gapIdx = 0; gapIdx < _gapDivisions; ++gapIdx) {
        // Y-position of this filament relative to gap center
        double filY;
        if (_gapDivisions == 1) {
            filY = 0.0;
        } else {
            filY = -(gapLength - gapGrid) / 2.0 + gapIdx * gapGrid;
        }

        // Absolute position of the "cross" current (into page)
        // For center-leg gap (gapX~0): place at the centerpost edge
        double crossX, crossY;
        if (std::abs(gapX) < 1e-10) {
            crossX = 0.0;
            crossY = gapY + filY;
        } else {
            crossX = gapX;
            crossY = gapY + filY;
        }

        // Mirror image "dot" current (out of page): reflected about x=0
        double dotX = -crossX;
        double dotY = crossY;

        // ---- Sum over all image units ----
        // This is the core of the method of images from Sullivan's code:
        // center_matrix = ones(size(y_temp))' * x_temp + j*y_temp' * ones(size(x_temp))
        // where x_temp = 2*hw * linspace(-unit_of_X, unit_of_X, ...)
        //       y_temp = bw  * linspace(-unit_of_Y, unit_of_Y, ...)
        for (int nx = -_imageUnitsX; nx <= _imageUnitsX; ++nx) {
            for (int ny = -_imageUnitsY; ny <= _imageUnitsY; ++ny) {
                double unitCenterX = nx * 2.0 * hw;
                double unitCenterY = ny * bw;

                // "Cross" current position in this image unit
                double srcCrossX = crossX + unitCenterX;
                double srcCrossY = crossY + unitCenterY;

                // "Dot" current position in this image unit
                double srcDotX = dotX + unitCenterX;
                double srcDotY = dotY + unitCenterY;

                // Biot-Savart for "cross" current (INTO page, +z direction)
                // B_x = +(mu_0*I)/(2*pi) * dy/r^2
                // B_y = -(mu_0*I)/(2*pi) * dx/r^2
                {
                    double dx = xP - srcCrossX;
                    double dy = yP - srcCrossY;
                    double r2 = dx * dx + dy * dy;
                    // Avoid division by zero (same approach as MATLAB code:
                    // Rp_abs = ((Rp_abs == 0) + Rp_abs); )
                    if (r2 < 1e-30) r2 = 1.0;

                    double coeff = u0 * I_per_div / (2.0 * std::numbers::pi * r2);
                    Bx_total += coeff * dy;
                    By_total -= coeff * dx;
                }

                // Biot-Savart for "dot" current (OUT OF page, -z direction)
                // Opposite sign current
                {
                    double dx = xP - srcDotX;
                    double dy = yP - srcDotY;
                    double r2 = dx * dx + dy * dy;
                    if (r2 < 1e-30) r2 = 1.0;

                    double coeff = u0 * (-I_per_div) / (2.0 * std::numbers::pi * r2);
                    Bx_total += coeff * dy;
                    By_total -= coeff * dx;
                }
            }
        }
    }

    // Convert B to H: H = B / mu_0. SIGN (ABT #1463, see the Roshen model): the field must be that
    // of a positive magnetizing MMF, whose flux runs -y down a lateral leg. As summed above, the
    // lateral gap's "cross" (at the leg, x = gapX) and mirrored "dot" (x = -gapX) filaments both
    // give +y in the window between them for H_g > 0, the opposite; hence the minus. A central
    // gap's cross and dot filaments both sit at x = 0 and cancel exactly, so this model returns
    // no field for a central-leg gap (a separate defect, not a sign).
    double Hx = -attenuationFactor * Bx_total / u0;
    double Hy = -attenuationFactor * By_total / u0;

    if (std::isnan(Hx) || std::isnan(Hy)) {
        throw NaNResultException("NaN found in Sullivan's fringing field model");
    }

    ComplexFieldPoint complexFieldPoint;
    complexFieldPoint.set_real(Hx);
    complexFieldPoint.set_imaginary(Hy);
    complexFieldPoint.set_point(inducedFieldPoint.get_point());
    if (inducedFieldPoint.get_turn_index()) {
        complexFieldPoint.set_turn_index(inducedFieldPoint.get_turn_index().value());
    }
    return complexFieldPoint;
}

// ============================================================================
// MagneticFieldStrengthImagedMmfSheetsModel (ABT #1409)
// ============================================================================

ComplexFieldPoint MagneticFieldStrengthImagedMmfSheetsModel::get_magnetic_field_strength_between_two_points(
    [[maybe_unused]] const FieldPoint& inducingFieldPoint,
    [[maybe_unused]] const FieldPoint& inducedFieldPoint,
    [[maybe_unused]] std::optional<size_t> inducingWireIndex) {
    throw CalculationException(ErrorCode::CALCULATION_ERROR,
        "IMAGED_MMF_SHEETS computes the window field of all turns and gaps together; it has no point-pair kernel");
}

// Uniform current `current` over [x1,x2] x [y1,y2], out of the plane. Exact 2D integral of the
// line-current field: with u = x - x', v = y - y' and P(u, v) = u ln(u^2 + v^2) + 2 v atan(u / v),
//   Hx = -J/(4 pi) [P]_{u,v},   Hy = J/(4 pi) [P(v, u)]_{u,v}   (double differences over both limits).
std::pair<double, double> MagneticFieldStrengthImagedMmfSheetsModel::rectangle_field(double x, double y, double x1, double x2, double y1, double y2, double current) {
    auto P = [](double u, double v) {
        double r2 = u * u + v * v;
        double logarithmicPart = r2 > 0 ? u * std::log(r2) : 0.0;
        double angularPart = v != 0 ? 2 * v * std::atan(u / v) : 0.0;
        return logarithmicPart + angularPart;
    };
    double density = current / ((x2 - x1) * (y2 - y1));
    double ua = x - x2, ub = x - x1, va = y - y2, vb = y - y1;
    double Hx = -density / (4 * std::numbers::pi) * (P(ub, vb) - P(ub, va) - P(ua, vb) + P(ua, va));
    double Hy = density / (4 * std::numbers::pi) * (P(vb, ub) - P(va, ub) - P(vb, ua) + P(va, ua));
    return {Hx, Hy};
}

// Uniform current `current` spread over the segment x = xs, y in [y1, y2].
std::pair<double, double> MagneticFieldStrengthImagedMmfSheetsModel::vertical_sheet_field(double x, double y, double xs, double y1, double y2, double current) {
    double u = x - xs;
    double va = y - y2, vb = y - y1;
    double linearDensity = current / (y2 - y1);
    double Hx = -linearDensity / (4 * std::numbers::pi) * (std::log(u * u + vb * vb) - std::log(u * u + va * va));
    double Hy;
    if (u != 0) {
        Hy = linearDensity / (2 * std::numbers::pi) * (std::atan(vb / u) - std::atan(va / u));
    }
    else if (va * vb > 0) {
        Hy = 0;
    }
    else {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "IMAGED_MMF_SHEETS: field requested on a gap sheet");
    }
    return {Hx, Hy};
}

std::pair<double, double> MagneticFieldStrengthImagedMmfSheetsModel::filament_field(double x, double y, double xs, double ys, double current) {
    double dx = x - xs, dy = y - ys;
    double r2 = dx * dx + dy * dy;
    if (r2 == 0) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "IMAGED_MMF_SHEETS: field requested on a filament");
    }
    return {-current * dy / (2 * std::numbers::pi * r2), current * dx / (2 * std::numbers::pi * r2)};
}

std::pair<double, double> MagneticFieldStrengthImagedMmfSheetsModel::source_field(const Source& source, double x, double y, double current) {
    if (source.x1 == source.x2 && source.y1 == source.y2) {
        return filament_field(x, y, source.x1, source.y1, current);
    }
    if (source.x1 == source.x2) {
        return vertical_sheet_field(x, y, source.x1, source.y1, source.y2, current);
    }
    return rectangle_field(x, y, source.x1, source.x2, source.y1, source.y2, current);
}

void MagneticFieldStrengthImagedMmfSheetsModel::setup(Magnetic magnetic, const std::vector<Wire>& wirePerWinding,
                                                      const std::vector<double>& currentDividerPerTurn,
                                                      const std::vector<int8_t>& currentDirectionPerWinding,
                                                      double frequency) {
    _turnSources.clear();
    _gapSources.clear();
    auto core = magnetic.get_core();
    if (core.get_shape_family() == CoreShapeFamily::T) {
        throw NotImplementedException("IMAGED_MMF_SHEETS: toroidal cores have no rectangular image frame");
    }
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException("IMAGED_MMF_SHEETS: core is not processed");
    }
    auto windingWindows = core.get_processed_description()->get_winding_windows();
    if (windingWindows.size() != 1) {
        throw NotImplementedException("IMAGED_MMF_SHEETS: only single-window cores are supported (this core has " +
                                      std::to_string(windingWindows.size()) + " winding windows)");
    }
    if (!windingWindows[0].get_width() || !windingWindows[0].get_height()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "IMAGED_MMF_SHEETS: the winding window has no width or height");
    }
    const double A = windingWindows[0].get_width().value();
    const double B = windingWindows[0].get_height().value();
    // Same frame as CoilMesherCenterModel for a single window: bounded left by the main column.
    const double frameLeftX = core.get_columns()[0].get_width() / 2;
    const double frameBottomY = -B / 2;

    int mirroringDimension = settings.get_magnetic_field_mirroring_dimension();
    if (mirroringDimension < 1) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "IMAGED_MMF_SHEETS needs the core walls (magnetic field mirroring dimension >= 1); without images the gap sheets are unbalanced");
    }
    bool completeCells = settings.get_magnetic_field_mirroring_complete_cells();
    double corePermeability = core.get_initial_permeability(defaults.ambientTemperature);
    int cellsX = mirroringDimension;
    int cellsY = mirroringDimension;
    if (completeCells) {
        cellsX = mirroringDimension * static_cast<int>(std::max(1.0, std::ceil(B / A)));
        cellsY = mirroringDimension * static_cast<int>(std::max(1.0, std::ceil(A / B)));
    }
    int mMinimum = completeCells ? -2 * cellsX - 1 : -cellsX;
    int mMaximum = completeCells ? 2 * cellsX : cellsX;
    int nMinimum = completeCells ? -2 * cellsY - 1 : -cellsY;
    int nMaximum = completeCells ? 2 * cellsY : cellsY;

    // The lattice rule and weights of CoilMesherCenterModel::generate_mesh_inducing_turn.
    auto appendImaged = [&](Source base, std::vector<Source>& out) {
        for (int m = mMinimum; m <= mMaximum; ++m) {
            for (int n = nMinimum; n <= nMaximum; ++n) {
                double multiplier = completeCells
                    ? std::pow((corePermeability - 1.0) / (corePermeability + 1.0), std::abs(m) + std::abs(n))
                    : (corePermeability - std::max(std::abs(m), std::abs(n))) / (corePermeability + std::max(std::abs(m), std::abs(n)));
                auto mapX = [&](double x) {
                    double a = x - frameLeftX;
                    return frameLeftX + ((m % 2 == 0) ? m * A + a : m * A + A - a);
                };
                auto mapY = [&](double y) {
                    double b = y - frameBottomY;
                    return frameBottomY + ((n % 2 == 0) ? n * B + b : n * B + B - b);
                };
                Source image = base;
                double xa = mapX(base.x1), xb = mapX(base.x2);
                double ya = mapY(base.y1), yb = mapY(base.y2);
                image.x1 = std::min(xa, xb);
                image.x2 = std::max(xa, xb);
                image.y1 = std::min(ya, yb);
                image.y2 = std::max(ya, yb);
                image.weight = base.weight * multiplier;
                image.isImage = (m != 0 || n != 0);
                out.push_back(image);
            }
        }
    };

    auto coil = magnetic.get_coil();
    if (!coil.get_turns_description()) {
        throw CoilNotProcessedException("IMAGED_MMF_SHEETS: the coil has no turns description");
    }
    const auto turns = coil.get_turns_description().value();
    if (currentDividerPerTurn.size() != turns.size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "IMAGED_MMF_SHEETS: " + std::to_string(currentDividerPerTurn.size()) +
                                    " current dividers for " + std::to_string(turns.size()) + " turns");
    }
    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        const auto& turn = turns[turnIndex];
        if (turn.get_additional_coordinates()) {
            throw NotImplementedException("IMAGED_MMF_SHEETS: turns with a second crossing (multi-column windings) are not supported");
        }
        size_t windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        auto wire = wirePerWinding.at(windingIndex);  // the dimension getters are non-const
        Source source;
        source.turnIndex = turnIndex;
        source.windingIndex = windingIndex;
        source.isImage = false;
        source.weight = currentDividerPerTurn[turnIndex] * currentDirectionPerWinding.at(windingIndex);
        double x = turn.get_coordinates()[0];
        double y = turn.get_coordinates()[1];
        if (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ) {
            source.x1 = source.x2 = x;
            source.y1 = source.y2 = y;
        }
        else {
            double halfWidth = wire.get_maximum_conducting_width() / 2;
            double halfHeight = wire.get_maximum_conducting_height() / 2;
            source.x1 = x - halfWidth;
            source.x2 = x + halfWidth;
            source.y1 = y - halfHeight;
            source.y2 = y + halfHeight;
        }
        appendImaged(source, _turnSources);
    }

    // Gap sheets: MMF per ampere of magnetizing current.
    auto gapping = core.get_functional_description().get_gapping();
    bool anyFunctional = false;
    for (const auto& gap : gapping) {
        if (gap.get_type() == GapType::SUBTRACTIVE || gap.get_type() == GapType::ADDITIVE) {
            anyFunctional = true;
        }
    }
    if (!anyFunctional) {
        return;
    }
    auto reluctanceModel = ReluctanceModel::factory();
    InitialPermeability initialPermeabilityModel;
    double initialPermeability = initialPermeabilityModel.get_initial_permeability(core.resolve_material(), std::nullopt, std::nullopt, frequency);
    double totalReluctance = reluctanceModel->get_core_reluctance(core, initialPermeability).get_core_reluctance();
    double numberTurns = static_cast<double>(coil.get_functional_description()[0].get_number_turns());
    double fluxPerAmpere = numberTurns / totalReluctance;

    auto halfSectionOf = [](const CoreGap& gap) {
        if (!gap.get_section_dimensions() || !gap.get_coordinates()) {
            throw GapException("IMAGED_MMF_SHEETS: a gap has no section dimensions or coordinates");
        }
        return gap.get_section_dimensions().value()[0] / 2;
    };
    std::map<long, double> lateralColumnReluctance;
    for (const auto& gap : gapping) {
        double x = gap.get_coordinates().value()[0];
        if (std::abs(x) >= halfSectionOf(gap)) {
            lateralColumnReluctance[std::lround(x * 1e6)] += reluctanceModel->get_gap_reluctance(gap).get_reluctance();
        }
    }
    double lateralConductance = 0;
    for (const auto& [column, reluctance] : lateralColumnReluctance) {
        lateralConductance += 1.0 / reluctance;
    }
    for (const auto& gap : gapping) {
        if (gap.get_type() != GapType::SUBTRACTIVE && gap.get_type() != GapType::ADDITIVE) {
            continue;  // residual mating surfaces carry no sheet (ABT #832)
        }
        double x = gap.get_coordinates().value()[0];
        if (x < 0) {
            continue;  // the other side's window
        }
        double halfSection = halfSectionOf(gap);
        double columnFluxPerAmpere = fluxPerAmpere;
        double sheetX;
        if (std::abs(x) < halfSection) {
            sheetX = x + halfSection;
        }
        else {
            columnFluxPerAmpere = fluxPerAmpere * (1.0 / lateralColumnReluctance.at(std::lround(x * 1e6))) / lateralConductance;
            sheetX = x - halfSection;
        }
        double mmfPerAmpere = columnFluxPerAmpere * reluctanceModel->get_gap_reluctance(gap).get_reluctance();
        double y = gap.get_coordinates().value()[1];
        Source sheet;
        sheet.x1 = sheet.x2 = sheetX;
        sheet.y1 = y - gap.get_length() / 2;
        sheet.y2 = y + gap.get_length() / 2;
        sheet.weight = -currentDirectionPerWinding.at(0) * mmfPerAmpere;
        sheet.windingIndex = 0;
        sheet.isImage = false;
        appendImaged(sheet, _gapSources);
    }
}

std::pair<double, double> MagneticFieldStrengthImagedMmfSheetsModel::turns_field(double x, double y, const std::vector<double>& currentPerWinding,
                                                                                 std::optional<size_t> excludedTurn) const {
    double Hx = 0, Hy = 0;
    for (const auto& source : _turnSources) {
        if (!source.isImage && excludedTurn && source.turnIndex == excludedTurn) {
            continue;
        }
        double current = source.weight * currentPerWinding.at(source.windingIndex);
        if (current == 0) {
            continue;
        }
        auto [hx, hy] = source_field(source, x, y, current);
        Hx += hx;
        Hy += hy;
    }
    return {Hx, Hy};
}

std::pair<double, double> MagneticFieldStrengthImagedMmfSheetsModel::gaps_field(double x, double y) const {
    double Hx = 0, Hy = 0;
    for (const auto& source : _gapSources) {
        auto [hx, hy] = source_field(source, x, y, source.weight);
        Hx += hx;
        Hy += hy;
    }
    return {Hx, Hy};
}

} // namespace OpenMagnetics
