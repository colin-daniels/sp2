#ifndef SP2_ENUMS_HPP
#define SP2_ENUMS_HPP

#include <string>
#include <unordered_map>
#include <src/common/io/util.hpp>

#ifdef SP2_CAN_HASH_ENUM
#define SPECIALIZE_ENUM_HASH(E)
#else
#define SPECIALIZE_ENUM_HASH(E)                                                \
namespace std {                                                                \
template<>                                                                     \
struct hash<E>                                                                 \
{                                                                              \
    typedef E argument_type;                                                   \
    typedef std::size_t result_type;                                           \
    result_type operator()(argument_type const& s) const                       \
    {                                                                          \
        using ut = std::underlying_type_t<E>;                                  \
        return std::hash<ut>{}(static_cast<ut>(s));                            \
    }                                                                          \
};                                                                             \
} // namespace std
#endif // SP2_CAN_HASH_ENUM

namespace sp2 {

template<class E>
using enum_map_t = std::pair<E, const char*>[];

/// variable template used for converting enums to/from strings, each new enum
/// that needs string conversion should explicitly specialize this
template<class E>
constexpr enum_map_t<E> enum_map = {};
// usage:
// enum class color {RED, BLUE};
//
// template<>
// constexpr enum_map_t<color>
//      enum_map<color> = {{color::RED, "red"}, {color::BLUE, "blue"}};

inline namespace {
/// convert enumeration E to string using enum_map<E>
template<class E>
std::string enum_to_str(E val);

/// convert std::string to enumeration E using enum_map<E> and default value
template<class E, bool case_sensitive = true>
E enum_from_str(const std::string &str, E default_value = static_cast<E>(0));

////////////////////////////////////////////////////////////////////////////////
// Implementations                                                            //
////////////////////////////////////////////////////////////////////////////////

template<class E>
std::string enum_to_str(E val)
{
    const static auto map = std::unordered_map<E, std::string>(
        std::begin(enum_map<E>), std::end(enum_map<E>));

    return map.at(val);
}


template<class E, bool case_sensitive>
E enum_from_str(const std::string &str, E default_value)
{
    const static auto map = []{
        auto result = std::unordered_map<std::string, E>{};
        for (auto v : enum_map<E>)
            result.emplace(io::tolower(v.second), v.first);

        return result;
    }();

    const auto loc = map.find(io::tolower(str));
    // return a default value if we can't find the string in the map or if
    // we CAN but the case-sensitive comparison doesn't check out
    bool not_found = loc == map.end() ||
        (case_sensitive && enum_to_str(loc->second) != str);

    return not_found ? default_value : loc->second;
}

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////
// Basic enum types used throughout the code                                  //
////////////////////////////////////////////////////////////////////////////////

enum class atom_type_old : int
{
    HYDROGEN = 0,
    CARBON   = 1
};

template<>
constexpr enum_map_t<atom_type_old> enum_map<atom_type_old> = {
    {atom_type_old::CARBON,   "C"},
    {atom_type_old::HYDROGEN, "H"}
};

enum class bond_type : int
{
    HYDROGEN_HYDROGEN = 0,
    HYDROGEN_CARBON   = 1,
    CARBON_CARBON     = 3
};

/// primary program task
enum class run_type : int
{
    NONE    = 0,
    RELAX   = 1,
    ATAC    = 2,
    SYMM    = 3,
    PHONOPY = 4
};

template<>
constexpr enum_map_t<run_type> enum_map<run_type> = {
    {run_type::RELAX,   "relax"},
    {run_type::ATAC,    "atac"},
    {run_type::SYMM,    "symm"},
    {run_type::PHONOPY, "phonopy"}
};

/// class of potential to use for force/energy calculation
enum class potential_type : int
{
    NONE   = 0,
    REBO   = 1,
    LAMMPS = 2,
    PHOSPHORENE = 3
};

template<>
constexpr enum_map_t<potential_type> enum_map<potential_type> = {
    {potential_type::REBO,   "rebo"},
    {potential_type::LAMMPS, "lammps"}
};

/// atom types for electrodes/transport
enum class el_type : int
{
    DEVICE        = 0,
    CONTACT_1_PL1 = 1,
    CONTACT_1_PL2 = 2,
    CONTACT_2_PL1 = 3,
    CONTACT_2_PL2 = 4
};

/// get bond type from atom types
static inline bond_type get_bond_type(const atom_type_old a, const atom_type_old b)
{
    return static_cast<bond_type>(
        (static_cast<int>(a) + 1)
        * (static_cast<int>(b) + 1) - 1
    );
}

} // namespace sp2

#include "atom_types.hpp"

SPECIALIZE_ENUM_HASH(sp2::atom_type)
SPECIALIZE_ENUM_HASH(sp2::run_type)
SPECIALIZE_ENUM_HASH(sp2::potential_type)

#endif // SP2_ENUMS_HPP
