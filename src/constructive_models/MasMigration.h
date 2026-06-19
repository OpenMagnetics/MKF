#pragma once

// Backwards-compatibility helper for loading pre-1.0 MAS documents.
//
// MAS 1.0 renamed many enum values: PD1/PD2/PD3 (was P1/P2/P3),
// I/II/III/IV (was OVC-I/...), camelCase for Application/ConnectionType/
// Topologies/FlybackModes/Market/etc. The auto-generated MAS.hpp
// from_json functions only accept the new spellings.
//
// This helper walks a json value in-place and rewrites old enum strings
// to their 1.0 equivalents. MKF wrapper-class constructors call it as
// the first line, so any pre-1.0 document loaded through those wrappers
// is transparently migrated before deserialization.
//
// The mapping is the C++ port of MAS/scripts/migrate-to-1.0.py and is
// the single source of truth on the C++ side. If you update the mapping
// in one place, update the other.

#include <nlohmann/json.hpp>

// ----------------------------------------------------------------------------
// Topology operating-point aliases (RFC 0006)
//
// Six topologies (forward, pushPull, llcResonant, isolatedBuck,
// isolatedBuckBoost, flybuck) inherit baseOperatingPoint via allOf with no
// extras. Quicktype merges them away into the base class (TopologyExcitation),
// so the *OperatingPoint names from earlier schema versions disappear.
// These aliases keep MKF source code that references them by name compiling.
//
// ----------------------------------------------------------------------------
#if __has_include(<MAS.hpp>)
#include <MAS.hpp>
namespace MAS {
    // RFC 0006: six topologies merged into TopologyExcitation base
    using BaseOperatingPoint              = TopologyExcitation;
    using ForwardOperatingPoint           = TopologyExcitation;
    using PushPullOperatingPoint          = TopologyExcitation;
    using LlcOperatingPoint               = TopologyExcitation;
    using IsolatedBuckOperatingPoint      = TopologyExcitation;
    using IsolatedBuckBoostOperatingPoint = TopologyExcitation;
    using FlybuckOperatingPoint           = TopologyExcitation;

    // PEAS integration (ae4414f): several types moved from MAS to PEAS utils.json
    // or renamed to avoid collisions with the broader PEAS vocabulary.
    using Application                  = MagneticApplication;
    using Topologies                   = Topology;
    using InsulationType               = IsolationClass;
    using InsulationCoordinationOutput = InsulationCoordination;
    using DistributorInfo              = Utils;
    // LossElementPerHarmonic (from PEAS outputs schema) replaces WindingLossElement
    using WindingLossElement           = LossElementPerHarmonic;

    // SubApplication was a MAS enum class; it is now a plain string field in
    // PEAS designRequirementsBase.  Define it as a namespace of string constants
    // so existing comparisons (value() == SubApplication::COMMON_MODE_NOISE_FILTERING)
    // continue to compile as std::string comparisons.
    namespace SubApplication {
        constexpr const char* COMMON_MODE_NOISE_FILTERING      = "commonModeNoiseFiltering";
        constexpr const char* DIFFERENTIAL_MODE_NOISE_FILTERING = "differentialModeNoiseFiltering";
        constexpr const char* POWER_FILTERING                   = "powerFiltering";
        constexpr const char* TRANSFORMING                      = "transforming";
        constexpr const char* ISOLATION                         = "isolation";
    }

}

#endif

// Convert the DesignRequirements.application string (PEAS-sourced, camelCase JSON value)
// to a MagneticApplication enum.  Throws if the string is not a valid application.
inline MAS::MagneticApplication magnetic_application_from_string(const std::string& s) {
    MAS::MagneticApplication x;
    nlohmann::json j = s;
    MAS::from_json(j, x);
    return x;
}

namespace OpenMagnetics::compat {

// Walk j recursively. Any string value that matches a pre-1.0 enum
// spelling is replaced by its 1.0 equivalent. Object keys are NOT
// rewritten (only values).
//
// Idempotent: running on a 1.0 document is a no-op.
void migrate_pre_1_0(nlohmann::json& j);

}  // namespace OpenMagnetics::compat
