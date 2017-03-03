#include "symm/system_control_t.hpp"
#include "symm/util.hpp"
#include "common/math/vec3_util.hpp"

#include <queue>

using namespace std;
using namespace sp2;

symm::system_control_t::system_control_t()
{
    unit_cell = 2;
    for (int i = 0; i < 3; ++i)
        lattice[i][i] = unit_cell;

    sub_sys.init(lattice, {}, {});

    // setup bond control for carbon
    bond_control.init(lattice, 2.0, 0.1);
}

void symm::system_control_t::set_group(const space_group_t &input) {
    group = input;}

structure_t symm::system_control_t::get_structure() const
{
    auto output = structure_t(lattice, types, v3tod(position));
    output.space_group = group.get_name();
    output.n_symm = group.n_symm();
    return output;
}

structure_t symm::system_control_t::get_full_structure() const
{
    auto output = structure_t(lattice, get_full_types(), v3tod(get_full_pos()));
    output.space_group = group.get_name();
    output.n_symm = group.n_symm();
    return output;
}

void symm::system_control_t::set_structure(const structure_t &input) {
    if (input.lattice[0][0] != unit_cell)
    {
        unit_cell = input.lattice[0][0];
        for (int i = 0; i < 3; ++i)
            lattice[i][i] = unit_cell;

        bond_control.set_lattice(lattice);
    }

    types = input.types;
    position = dtov3(input.positions);
}

std::vector<double> symm::system_control_t::get_gradient() const {
    return v3tod(gradient);}

double symm::system_control_t::get_value() const {
    return total_potential;}

/// calculate forces and total potential for the current structure
void symm::system_control_t::update()
{
    constexpr bool limited = true;

    // input new structure into the rebo system to calculate forces/potential
    structure_t structure = limited ?
                            get_limited_structure() : get_full_structure();

    sub_sys.set_structure(structure);

    // calculate potential/force
    sub_sys.update();

    // get the total potential
    total_potential = 0;
    if (limited)
    {
        auto potentials = sub_sys.get_potentials();
        for (size_t i = 0; i < position.size(); ++i)
            total_potential += 0.5 * potentials[i];
    }
    else
    {
        total_potential = sub_sys.get_value() / group.n_symm();
    }

    // get the gradient
    gradient = dtov3(sub_sys.get_gradient());
    gradient.resize(position.size());

    // cout << total_potential << endl;

    // gyroid surface penalty potential
    for (size_t i = 0; i < position.size(); ++i)
    {
        constexpr double scale = 0.7;

        vec3_t deriv = {},
            pos = position[i] / unit_cell;

        total_potential += scale * gyroid_ptnl(pos, deriv);
        gradient[i] += scale * deriv;
    }

//     // get unit cell force
//     const auto &bc = sub_sys.get_bond_control();
//
//     auto graph = bc.get_graph();
//     auto deltas = bc.get_bond_deltas(),
//         forces = sub_sys.get_bond_forces();
//
//     force.resize(na_base * 3 + 1, 0);
//
//     uc_force = -blas_ddot(deltas.size(), deltas.data(), 1, bond_force.data(), 1) / (2 * unit_cell);
//     if (unit_cell == uc_min && uc_force < 0)
//         uc_force = 0;
//     else if (unit_cell == uc_max && uc_force > 0)
//         uc_force = 0;
// //    force.back() = uc_force;
}

structure_t symm::system_control_t::get_limited_structure(size_t max_dist)
{
    auto full_pos = get_full_pos();
    auto full_types = get_full_types();

    // get the bond graph for the structure
    bond_control.update(v3tod(full_pos));

    auto offsets = bond_control.get_bond_offsets(),
        bond_ids = bond_control.get_bond_ids();

    auto graph = bond_control.get_graph();

    // breadth first search to select atoms withing max_dist bonds of the
    // irreducible representation
    queue<int> to_visit;
    vector<bool> visited(full_pos.size(), false);

    for (size_t i = 0; i < position.size(); ++i)
    {
        to_visit.push(i);
        visited[i] = true;
    }

    for (size_t i = 0; i + 1 < max_dist; ++i)
    {
        size_t cur_size = to_visit.size();
        for (size_t j = 0; j < cur_size; ++j)
        {
            auto id = to_visit.front();
            to_visit.pop();

            for (auto neigh_id : graph.neighbors(id))
            {
                if (!visited[neigh_id])
                {
                    visited[neigh_id] = true;
                    to_visit.push(neigh_id);
                }
            }
        }
    }

    // only keep atoms which have been marked
    vector<double> selected_pos;
    vector<atom_type> selected_types;

    selected_pos.reserve(full_pos.size() * 3);
    selected_types.reserve(full_types.size());

    for (size_t i = 0; i < visited.size(); ++i)
    {
        if (!visited[i])
            continue;

        selected_types.push_back(full_types[i]);
        for (auto d : full_pos[i])
            selected_pos.push_back(d);
    }

    return structure_t(lattice, selected_types, selected_pos);
}

std::vector<vec3_t> symm::system_control_t::get_full_pos() const
{
    // scale down
    vector<vec3_t> full_pos = position;
    for (auto &v : full_pos)
        v /= unit_cell;

    // apply the symmetries
    full_pos = group.apply_symm(full_pos);

    // scale back up
    for (auto &v : full_pos)
        v *= unit_cell;

    return full_pos;
}


std::vector<atom_type> symm::system_control_t::get_full_types() const
{
    vector<atom_type> full_types;
    for (size_t i = 0; i < group.n_symm(); ++i)
        full_types.insert(full_types.end(), types.begin(), types.end());

    return full_types;
}
