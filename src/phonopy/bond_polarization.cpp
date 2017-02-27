#include <common/graph/ud_graph_t.hpp>
#include <common/vec3_t.hpp>
#include <airebo/system_control_t.hpp>
#include "bond_polarization.hpp"

#include "common/mat3x3_t.hpp"

double sp2::phonopy::raman_intensity(
    double frequency, double temperature,
    vec3_t incident, vec3_t scattered,
    const std::vector<vec3_t> &eigs,
    const graph::ud_graph_t &bond_graph,
    const std::vector<vec3_t> &bonds,
    const std::vector<atom_types> &types,
    const std::unordered_map<bond_types, pol_constant_t> &pol_constants
)
{
    constexpr double hk = 0.228988; // hbar / k_b in [K][cm]

    double distribution = 0;
    if (temperature != 0)
        distribution = 1 / (std::exp(hk * frequency / temperature) - 1);

    mat3x3_t tensor = raman_tensor(
        eigs, bond_graph, bonds, types, pol_constants
    );
    const double sum = dot(incident, scattered.mul_3x3(tensor));

    return ((distribution + 1) / frequency) * (sum * sum);
}

sp2::mat3x3_t sp2::phonopy::raman_tensor(
    const std::vector<vec3_t> &eigs,
    const graph::ud_graph_t &bond_graph,
    const std::vector<vec3_t> &bonds,
    const std::vector<atom_types> &types,
    const std::unordered_map<bond_types, pol_constant_t> &pol_constants
)
{
    // kronecker delta value
    constexpr double kdelta[3][3] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };

    mat3x3_t tensor{};
    for (auto edge : bond_graph.edges())
    {
        const int bond_id = edge.id;

        // phonon eigenvector for this atom
        const vec3_t &eig = eigs[edge.a];

        // unit bond vector and length, used later
        const double len = bonds[bond_id].mag();
        const vec3_t rhat = bonds[bond_id] / len;

        const auto &pc = pol_constants.at(btype(types[edge.a], types[edge.b]));
        double const_one = pc.c1, // a || -   a |-
              dconst_one = pc.c2, // a'|| -   a'|-
              dconst_two = pc.c3; // a'|| + 2 a'|-

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                tensor[i][j] -= dot(rhat, eig) * (
                      (dconst_two / 3) * kdelta[i][j]
                    + dconst_one * (rhat[i] * rhat[j] - kdelta[i][j] / 3)
                ) + (const_one / len) * (
                      (rhat[i] * eig[j] + rhat[j] * eig[i])
                    - 2 * rhat[i] * rhat[j] * dot(rhat, eig)
                );
            }
        }
    }

    return tensor;
}

std::vector<std::pair<double, double>> sp2::phonopy::raman_spectra(
    vec3_t incident, vec3_t scattered, double temperature,
    const std::vector<std::pair<double, std::vector<vec3_t>>> &modes,
    const structure_t &structure
)
{
    sp2::airebo::system_control_t sys;
    sys.init(structure);
    sys.update();

    auto bond_graph = sys.get_bond_control().get_graph();
    auto bond_deltas = dtov3(sys.get_bond_control().get_bond_deltas());

    // convert from old to new types for now
    std::vector<atom_types> types;
    for (auto t : structure.types)
        types.push_back(t == atom_type::CARBON ? atype::C : atype::H);

    std::vector<std::pair<double, double>> result;
    for (auto mode : modes)
    {
        double frequency = mode.first;
        const auto &eigs = mode.second;

        // calculate raman intensity for given the incident/scattered light
        // polarization directions and the current eigenmode
        double intensity = raman_intensity(frequency, temperature,
            incident, scattered, eigs, bond_graph, bond_deltas, types);

        result.emplace_back(frequency, intensity);
    }

    // note that we don't normalize to unity, the user must
    return result;
}
