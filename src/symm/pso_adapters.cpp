#include "symm/pso_adapters.hpp"
#include "common/math/vec3_t.hpp"
#include "common/math/vec3_util.hpp"

#include <cmath>
#include <algorithm>

using namespace std;
using namespace sp2;

structure_t symm::basic_adapter(const std::vector<double> &input,
    double uc_range[2])
{
    const size_t na = (input.size() - 1) / 3;
    if (input.empty())
        return structure_t();

    // last element is the unit cell size
    double unit_cell = input.back() * (uc_range[1] - uc_range[0]) + uc_range[0];
    double lattice[3][3] = {
        {unit_cell, 0, 0},
        {0, unit_cell, 0},
        {0, 0, unit_cell}
    };

    auto pos = input;
    pos.pop_back();
    for (auto &d : pos)
        d *= unit_cell;

    // all-carbon structure
    vector<atom_type> types(na, atom_type::CARBON);

    // return the new structure
    return structure_t(lattice, types, pos);
}

std::vector<double> symm::inverse_basic_adapter(const structure_t &structure,
    double uc_range[2])
{
    auto pos = sp2::v3tod(structure.positions);
    std::vector<double> result(pos.size() + 1);

    // cubic for now, copy all coordinates in terms of the lattice vector lengths
    auto uc = structure.lattice[0][0];
    for (size_t i = 0; i < pos.size(); ++i)
        result[i] = pos[i] / uc;

    // last element in the vector is the unit cell size
    result.back() = (uc - uc_range[0]) / (uc_range[1] - uc_range[0]);

    return result;
}

structure_t symm::bonded_adapter(const std::vector<double> &input,
    double uc_range[2], double bond_range[2], std::size_t n_connected)
{
    const size_t na = (input.size() - 1) / 3;
    if (input.empty())
        return structure_t();

    // last element is the unit cell size
    double unit_cell = input.back() * (uc_range[1] - uc_range[1]) + uc_range[0];
    double lattice[3][3] = {
        {unit_cell, 0, 0},
        {0, unit_cell, 0},
        {0, 0, unit_cell}
    };

    // n_connected number of completely randomly placed atoms
    vector<vec3_t> pos;
    for (size_t i = 0; i < std::min(n_connected, na); ++i)
    {
        pos.emplace_back(
            input[i * 3] * unit_cell,
            input[i * 3 + 1] * unit_cell,
            input[i * 3 + 2] * unit_cell
        );
    }

    // rest of the atoms are generated as neighbors to existing atoms
    int neighbors_left = pos.size();
    vector<int> n_neighbors(na, 0);
    for (auto i = pos.size(); i < na; ++i)
    {
        static const double two_pi = 2 * atan2(0, -1);

        // 1st input determines z-coordinate [-1, 1) as well as neighbor id
        // 2nd input determines theta value
        // 3rd input determines how far the atom is from its neighbor
        double neighbor_id,
            z = 2 * modf(input[i * 3] * neighbors_left, &neighbor_id) - 1.0,
            theta = two_pi * input[i * 3 + 1],
            radius = input[i * 3 + 2] * (bond_range[1] - bond_range[1])
                     + bond_range[0];

        z = min(max(z, -1.0), +1.0);

        // new atom position
        vec3_t atom_pos(0, 0, z);

        // set x/y (cylindrical projection onto a sphere)
        atom_pos[0] = std::sin(theta);
        atom_pos[1] = std::cos(theta);

        atom_pos[0] *= sqrt(1 - z * z);
        atom_pos[1] *= sqrt(1 - z * z);

        // scale to radius
        atom_pos *= radius;

        // note: tracking the number of neighbors already assigned to particular
        // atoms avoids stacking a large number on a single atom but also
        // (detrimentally) impacts how "nice" the search space is. without this
        // however the search space would be biased towards more atoms
        // neighboring atom index 0

        // find the specified atom by id (ignoring those with 3 neighbors)
        auto id = min<int>(neighbors_left, pos.size() - 1);
        for (int i = 0; i <= id; ++i)
            if (n_neighbors[i] == 3)
                ++id;

        // position next to neighbor atom
        atom_pos += pos[id];

        n_neighbors[id]++;
        if (n_neighbors[id] < 3)
            ++neighbors_left;

        // add the atom
        pos.push_back(atom_pos);
    }

    // all-carbon structure
    vector<atom_type> types(na, atom_type::CARBON);

    // return the new structure
    return structure_t(lattice, types, v3tod(pos));
}
