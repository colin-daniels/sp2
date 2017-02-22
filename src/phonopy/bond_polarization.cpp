#include <common/graph/ud_graph_t.hpp>
#include <common/vec3_t.hpp>
#include "bond_polarization.hpp"

#include "common/mat3x3_t.hpp"

double sp2::phonopy::raman_intensity(
    vec3_t incident, vec3_t scattered,
    const std::vector<vec3_t> &eigs,
    const graph::ud_graph_t &bond_graph,
    const std::vector<vec3_t> &bonds,
    const std::vector<atom_types> &types,
    const std::unordered_map<bond_types, pol_constant_t> &pol_constants
)
{
    mat3x3_t tensor = raman_tensor(
        eigs, bond_graph, bonds, types, pol_constants
    );

    double sum = dot(incident, scattered.mul_3x3(tensor));
    return sum * sum;
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