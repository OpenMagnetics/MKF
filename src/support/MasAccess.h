#pragma once

// Read access, without a copy, to the std::optional members of the generated MAS classes.
//
// quicktype returns every std::optional member BY VALUE (`std::optional<Waveform> get_waveform() const`), so each
// presence test or read of an optional waveform, harmonics record, ... copies it (a 16k-sample Waveform is ~128 KB per
// call). MAS.hpp is generated and is never edited or post-processed; the members it declares private are reached here
// instead, through a pointer to the member that an explicit instantiation is allowed to name:
//
//     [temp.explicit]/12: "The usual access checking rules do not apply to names in a declaration of an explicit
//     instantiation." So `template struct Rob<Tag, &Class::member>;` may name a private member. The instantiation
//     defines a friend function that returns that pointer, and `object.*pointer` is an ordinary member access.
//
// Nothing here changes a class, a layout or a copy: it only reads. Only the (class, member) pairs MKF reads through a
// reference are declared; a pair that is not declared does not compile (there is no overload), so an unlisted access is
// a build error, never a silent copy. To read another member, add one MKF_MAS_OPTIONAL_REF line below.
//
// usage:  const auto& waveform = mas_access::waveform_ref(signalDescriptor);   // const std::optional<Waveform>&
//         if (mas_access::time_ref(waveform)) { ... }

#include <MAS.hpp>

#include <map>
#include <optional>
#include <string>
#include <vector>

// The macro must be expanded at namespace scope: the friend declaration it makes is found by argument-dependent
// lookup on the tag, and a friend function cannot be declared from a block scope.
#define MKF_MAS_OPTIONAL_REF(Class, member, ...)                                                                      \
    namespace mas_access_detail {                                                                                     \
    struct Class##_##member##_tag {                                                                                   \
        using pointer = __VA_ARGS__ MAS::Class::*;                                                                    \
        friend pointer mas_access_member(Class##_##member##_tag);                                                     \
    };                                                                                                                \
    template <Class##_##member##_tag::pointer Member>                                                                 \
    struct Class##_##member##_rob {                                                                                   \
        friend Class##_##member##_tag::pointer mas_access_member(Class##_##member##_tag) { return Member; }          \
    };                                                                                                                \
    template struct Class##_##member##_rob<&MAS::Class::member>;                                                      \
    }                                                                                                                 \
    inline const __VA_ARGS__& member##_ref(const MAS::Class& object) {                                                \
        return object.*mas_access_member(mas_access_detail::Class##_##member##_tag{});                                \
    }

namespace OpenMagnetics::mas_access {

MKF_MAS_OPTIONAL_REF(Waveform, time, std::optional<std::vector<double>>)

MKF_MAS_OPTIONAL_REF(SignalDescriptor, waveform, std::optional<MAS::Waveform>)
MKF_MAS_OPTIONAL_REF(SignalDescriptor, harmonics, std::optional<MAS::Harmonics>)
MKF_MAS_OPTIONAL_REF(SignalDescriptor, processed, std::optional<MAS::ProcessedWaveform>)

MKF_MAS_OPTIONAL_REF(OperatingPointExcitation, current, std::optional<MAS::SignalDescriptor>)
MKF_MAS_OPTIONAL_REF(OperatingPointExcitation, voltage, std::optional<MAS::SignalDescriptor>)
MKF_MAS_OPTIONAL_REF(OperatingPointExcitation, magnetizing_current, std::optional<MAS::SignalDescriptor>)
MKF_MAS_OPTIONAL_REF(OperatingPointExcitation, magnetic_flux_density, std::optional<MAS::SignalDescriptor>)

MKF_MAS_OPTIONAL_REF(PermeabilityPoint, frequency, std::optional<double>)
MKF_MAS_OPTIONAL_REF(PermeabilityPoint, temperature, std::optional<double>)
MKF_MAS_OPTIONAL_REF(PermeabilityPoint, tolerance, std::optional<double>)
MKF_MAS_OPTIONAL_REF(PermeabilityPoint, magnetic_field_dc_bias, std::optional<double>)
MKF_MAS_OPTIONAL_REF(PermeabilityPoint, magnetic_field_peak, std::optional<double>)
MKF_MAS_OPTIONAL_REF(PermeabilityPoint, magnetic_flux_density_peak, std::optional<double>)
MKF_MAS_OPTIONAL_REF(PermeabilityPoint, modifiers,
                     std::optional<std::map<std::string, MAS::InitialPermeabilitModifier>>)

MKF_MAS_OPTIONAL_REF(MagneticCore, processed_description, std::optional<MAS::CoreProcessedDescription>)

MKF_MAS_OPTIONAL_REF(Coil, turns_description, std::optional<std::vector<MAS::Turn>>)

} // namespace OpenMagnetics::mas_access

#undef MKF_MAS_OPTIONAL_REF
