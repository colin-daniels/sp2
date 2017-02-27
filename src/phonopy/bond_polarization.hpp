#ifndef SP2_BOND_POLARIZATION_HPP
#define SP2_BOND_POLARIZATION_HPP

#include "common/mat3x3_t.hpp"
#include "common/atom_types.hpp"
#include "common/structure_t.hpp"

#include <unordered_map>

namespace sp2 {
namespace phonopy {

struct pol_constant_t
{
    // default values are for C-C bonds
    double c1 = 0.32, // a || -   a |-
        c2 = 2.60,    // a'|| -   a'|-
        c3 = 7.55;    // a'|| + 2 a'|-
};

mat3x3_t raman_tensor(
    const std::vector<vec3_t> &eigs,
    const graph::ud_graph_t &bond_graph,
    const std::vector<vec3_t> &bonds,
    const std::vector<atom_types> &types,
    const std::unordered_map<bond_types, pol_constant_t> &pol_constants
);

double raman_intensity(
    double frequency, double temperature,
    vec3_t incident, vec3_t scattered,
    const std::vector<vec3_t> &eigs,
    const graph::ud_graph_t &bond_graph,
    const std::vector<vec3_t> &bonds,
    const std::vector<atom_types> &types,
    const std::unordered_map<bond_types, pol_constant_t>
        &pol_constants = {
            {btype(atype::C, atype::C), pol_constant_t{}},
            {btype(atype::C, atype::H), pol_constant_t{}}
        }
);

/// calculate raman spectra for a given system
/// \param incident polarization direction unit vector for the incident light
/// \param scattered polarization direction unit vector for the scattered light
/// \param modes vector of {frequency, mass-normalized eigenvectors} pairs
/// \param structure input structure
/// \return vector of pairs of {frequency [cm^-1], intensity [arb units]}
///         note: not normalized to maximum intensity
std::vector<std::pair<double, double>> raman_spectra(
    vec3_t incident, vec3_t scattered, double temperature,
    const std::vector<std::pair<double, std::vector<vec3_t>>> &modes,
    const structure_t &structure
);

} // namespace phonopy
} // namespace sp2

#endif // SP2_BOND_POLARIZATION_HPP
