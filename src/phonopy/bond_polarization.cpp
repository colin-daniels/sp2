#include <common/graph/ud_graph_t.hpp>
#include <common/vec3_t.hpp>
#include "bond_polarization.hpp"

namespace sp2 {

struct mat3x3_t
{
    double data[3][3];

    constexpr double const (&operator[](int i) const)[3]
    {
        return data[i];
    }

    constexpr double (&operator[](int i))[3]
    {
        return data[i];
    }

    operator const decltype(data)&()
    {
        return data;
    }
};

mat3x3_t raman_tensor(const std::vector<sp2::vec3_t> eigs,
    const sp2::graph::ud_graph_t &bond_graph,
    const std::vector<sp2::vec3_t> &bonds);

double raman_intensity(vec3_t incident, vec3_t scattered,
    const std::vector<sp2::vec3_t> eigs,
    const sp2::graph::ud_graph_t &bond_graph,
    const std::vector<sp2::vec3_t> &bonds
)
{
    mat3x3_t tensor = raman_tensor(eigs, bond_graph, bonds);
    double sum = dot(incident, scattered.mul_3x3(tensor));
    return sum * sum;
}


} // namespace sp2

sp2::mat3x3_t sp2::raman_tensor(const std::vector<sp2::vec3_t> eigs,
    const sp2::graph::ud_graph_t &bond_graph,
    const std::vector<sp2::vec3_t> &bonds)
{
    // kronecker delta value
    constexpr double kdelta[3][3] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };

    double const_one = 0.32, // a || -   a |-
          dconst_one = 2.60, // a'|| -   a'|-
          dconst_two = 7.55; // a'|| + 2 a'|-

    mat3x3_t tensor{};
    for (auto edge : bond_graph.edges())
    {
        const int atom_id = edge.a,
            bond_id = edge.id;

        // phonon eigenvector for this atom
        const vec3_t &eig = eigs[atom_id];

        // unit bond vector and length, used later
        const double len = bonds[bond_id].mag();
        const vec3_t rhat = bonds[bond_id] / len;

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