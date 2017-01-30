#ifndef SP2_ENUMS_HPP
#define SP2_ENUMS_HPP

#include <string>

namespace sp2 {
////////////////////////////////////////////////////////////////////////////////
// Basic enum types used throughout the code                                  //
////////////////////////////////////////////////////////////////////////////////

enum class atom_type : int
{
    HYDROGEN = 0,
    CARBON   = 1
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
    NONE     = 0,
    MINIMIZE = 1,
    ATAC     = 2,
    SYMM     = 3,
    PHONOPY  = 4
};

/// class of potential to use for force/energy calculation
enum class potential_type : int
{
    NONE   = 0,
    REBO   = 1,
    LAMMPS = 2
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

////////////////////////////////////////////////////////////////////////////////
// Enum conversions/misc                                                      //
////////////////////////////////////////////////////////////////////////////////

template<class T>
std::string enum_to_str(T val);

template<class T>
T enum_from_str(std::string str);

#ifndef BASIC_ENUMS_HPP_IMPL
extern template atom_type enum_from_str(std::string);
extern template std::string enum_to_str(atom_type);
#endif // BASIC_ENUMS_HPP_IMPL

/// get bond type from atom types
static inline bond_type get_bond_type(const atom_type a, const atom_type b)
{
    return static_cast<bond_type>(
        (static_cast<int>(a) + 1)
        * (static_cast<int>(b) + 1) - 1
    );
}

} // namespace sp2

#endif // SP2_ENUMS_HPP
