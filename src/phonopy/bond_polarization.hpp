#ifndef SP2_BOND_POLARIZATION_HPP
#define SP2_BOND_POLARIZATION_HPP

#include "common/math/mat3x3_t.hpp"
#include "common/atom_types.hpp"
#include "common/structure_t.hpp"

#include <unordered_map>

namespace sp2 {
namespace phonopy {

struct pol_constant_t
{
    double c1 = 0, ///< a || -   a |-
           c2 = 0, ///< a'|| -   a'|-
           c3 = 0; ///< a'|| + 2 a'|-

    /// maximum bond length
    double max_len = 0;

    pol_constant_t() = default;
};

template<bond_types>
constexpr pol_constant_t pol_const;

template<>
constexpr auto pol_const<bond_types::CC> =
    pol_constant_t{0.32, 2.60, 7.55, /* max_len */ 1.6};

template<>
constexpr auto pol_const<bond_types::CH> =
    pol_constant_t{0.32, 2.60, 7.55, /* max_len */ 1.3};

template<>
constexpr auto pol_const<bond_types::HH> =
//    pol_constant_t{0.32, 2.60, 7.55, /* max_len */ 2.1};
    pol_constant_t{0, 0, 0, /* max_len */ 1.1};

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
            {bond_types::CC, pol_const<bond_types::CC>},
            {bond_types::CH, pol_const<bond_types::CH>},
            {bond_types::HH, pol_const<bond_types::HH>}
        }
);

double raman_intensity_avg(
    bool backscatter,
    double frequency, double temperature,
    const std::vector<vec3_t> &eigs,
    const graph::ud_graph_t &bond_graph,
    const std::vector<vec3_t> &bonds,
    const std::vector<atom_types> &types,
    const std::unordered_map<bond_types, pol_constant_t>
        &pol_constants = {
            {bond_types::CC, pol_const<bond_types::CC>},
            {bond_types::CH, pol_const<bond_types::CH>},
            {bond_types::HH, pol_const<bond_types::HH>}
        }
);

/// calculate raman spectra for a given system
/// \param incident polarization direction unit vector for the incident light
/// \param scattered polarization direction unit vector for the scattered light
/// \param modes vector of {frequency, orthonormal mode eigenvectors} pairs
/// \param structure input structure
/// \return vector of pairs of {frequency [cm^-1], intensity [arb units]}
///         note: not normalized to maximum intensity
std::vector<std::pair<double, double>> raman_spectra(
    vec3_t incident, vec3_t scattered, double temperature,
    const std::vector<std::pair<double, std::vector<vec3_t>>> &modes,
    const structure_t &structure
);


std::vector<std::pair<double, double>> raman_spectra_avg(
    bool backscatter,
    double temperature,
    const std::vector<std::pair<double, std::vector<vec3_t>>> &modes,
    const structure_t &structure
);

} // namespace phonopy
} // namespace sp2

#endif // SP2_BOND_POLARIZATION_HPP
