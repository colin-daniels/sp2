#include <common/neighbor/bond_control_t.hpp>
#include "modeling.hpp"
#include "common/math/vec3_t.hpp"
#include "common/math/vec3_util.hpp"

sp2::structure_t sp2::util::construct_supercell(const sp2::structure_t &input,
    int supercell_dim[3])
{
    int n_rep = supercell_dim[0]
                * supercell_dim[1]
                * supercell_dim[2];

    structure_t supercell = input;

    // lengthen lattice vectors
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            supercell.lattice[i][j] *= supercell_dim[i];

    // get repeated types
    supercell.types.reserve(n_rep * input.types.size());
    for (int i = 1; i < n_rep; ++i)
    {
        supercell.types.insert(
            supercell.types.end(),
            input.types.begin(),
            input.types.end()
        );
    }

    // repeat positions
    vec3_t lattice[3] = {
        vec3_t(input.lattice[0]),
        vec3_t(input.lattice[1]),
        vec3_t(input.lattice[2])
    };

    auto original_pos = sp2::dtov3(input.positions),
        new_pos = std::vector<vec3_t>();

    auto add_new = [&](int i, int j, int k) {
        auto to_add = original_pos;

        for (auto &v : to_add)
            v += i * lattice[0] +
                 j * lattice[1] +
                 k * lattice[2];

        new_pos.insert(
            new_pos.end(),
            to_add.begin(),
            to_add.end()
        );
    };

    for (int i = 0; i < supercell_dim[0]; ++i)
        for (int j = 0; j < supercell_dim[1]; ++j)
            for (int k = 0; k < supercell_dim[2]; ++k)
                add_new(i, j, k);


    supercell.positions = sp2::v3tod(new_pos);

    return supercell;
}

sp2::structure_t sp2::util::construct_graphene(int m, int n, bool x_aligned)
{


    return sp2::structure_t();
}

sp2::structure_t sp2::util::make_hydrogen_terminated(
    const sp2::structure_t &input)
{
    // maximum distance between two atoms to be considered bonded in Angstroms
    constexpr double max_dist = 1.6,
        ch_dist = 1.09; // initial CH bond length for added atoms

    // use fbc::bond_control_t to calculate all the bonds between atoms
    // in the system
    fbc::bond_control_t bond_control;
    bond_control.init(input.lattice, max_dist, 0);
    bond_control.update(input.positions);

    auto types = input.types;
    auto pos = sp2::dtov3(input.positions);

    auto graph = bond_control.get_graph();
    auto bond_deltas = sp2::dtov3(bond_control.get_bond_deltas());

    for (auto idx : graph.vertices())
    {
        auto type = types[idx];
        if (type != atom_type::CARBON)
            continue;

        auto n_bonds = graph.degree(idx);
        if (n_bonds >= 3)
            continue;

        vec3_t bonds[3];
        int neigh_ids[3];

//        for (int j = 0; j < n_bonds; ++j)
//        {
//            auto edge = graph.edge(idx, j);
//            bonds[j] = bond_deltas[edge.id];
//            neigh_ids[j] = edge.b;
//        }
//
//        switch (n_bonds)
//        {
//        case 0:
//        case 1:
//        case 2:
//            break;
//        }
    }

//        int num_total = -1;
//        int per_atom_max = -1;
//
//        // update bond values
//        update();
//
//        // get references that will be needed later
//        const vector<double> &lengths = bond_control.get_bond_lengths(),
//            &deltas = bond_control.get_bond_deltas();
//        const vector<int> &bond_ids = bond_control.get_bond_ids(),
//            &offsets = bond_control.get_bond_offsets();
//
//        // rank atoms by number of bonds
//        vector<pair<double, int> > cvec;
//        for (int i = 0; i < na; ++i)
//            cvec.push_back(pair<double, int>(Na_vars[i * 4], i));
//        sort(cvec.begin(), cvec.end());
//
//        // add hydrogen
//        int h_left = num_total;
//
//        for (int i = 0; i < na; ++i)
//        {
//            int id = cvec[i].second;
//            if (types[id] == atom_type::HYDROGEN)
//                continue;
//
//            if (Na_vars[id * 4] >= 3 && h_left <= 0)
//                break;
//
//            int c_hyd = 0;
//            vector<pair<double, int> > bonds;
//            for (int j = offsets[id]; j < offsets[id + 1]; ++j)
//            {
//                bonds.push_back(pair<double, int>(cutoff[j], j));
//                if (types[bond_ids[j]] == atom_type::HYDROGEN && cutoff[j] >= 1)
//                    ++c_hyd;
//            }
//
//            if (per_atom_max > 0 && c_hyd >= per_atom_max)
//                continue;
//
//            int nbonds = min((int)bonds.size(), 2);
//            sort(bonds.begin(), bonds.end());
//
//            vector<int> added_ids;
//            for (int j = 0; j < nbonds; ++j)
//                added_ids.push_back(bonds[j].second);
//
//            // add hydrogen
//            if (nbonds == 2)
//            {
//                int b_id_a = added_ids[0],
//                    b_id_b = added_ids[1];
//
//                double avg[3], len = 0;
//                for (int j = 0; j < 3; ++j)
//                {
//                    avg[j] = deltas[b_id_a * 3 + j] / lengths[b_id_a]
//                             + deltas[b_id_b * 3 + j] / lengths[b_id_b];
//                    len += avg[j] * avg[j];
//                }
//
//                if (len == 0)
//                {
//                    auto avg_v = unit_normal(vec3_t(&deltas[b_id_a * 3]));
//                    for (int k = 0; k < 3; ++k)
//                        avg[k] = avg_v[k];
//
//                    len = 1;
//                }
//
//                for (int j = 0; j < 3; ++j)
//                    position.push_back(position[id * 3 + j] - 1.09 * avg[j]
//                                                              / sqrt(len));
//                types.push_back(atom_type::HYDROGEN);
//                --h_left;
//            }
//            else if (nbonds == 1)
//            {
//                vec3_t bond_ref(&deltas[added_ids[0] * 3]);
//                vec3_t axis = unit_normal(bond_ref);
//
//                double rot[3][3] = {};
//                const double pi = atan2(0, -1);
//
//                util::gen_rotation(axis, 2 * pi / 3, rot);
//                vec3_t bond_a = bond_ref.mul_3x3(rot).unit_vector() * 1.09;
//
//                util::gen_rotation(axis, 4 * pi / 3, rot);
//                vec3_t bond_b = bond_ref.mul_3x3(rot).unit_vector() * 1.09;
//
//                for (int j = 0; j < 3; ++j)
//                    position.push_back(position[id * 3 + j] + 1.09 * bond_a[j]);
//                types.push_back(atom_type::HYDROGEN);
//                --h_left;
//
//                if (per_atom_max > 1 && c_hyd == 0)
//                {
//                    for (int j = 0; j < 3; ++j)
//                        position.push_back(position[id * 3 + j] + 1.09 * bond_b[j]);
//                    types.push_back(atom_type::HYDROGEN);
//                    --h_left;
//                }
//            }
//
//            if (h_left == 0)
//                break;
//        }
//
//        // set number of bonds and atoms
//        na = types.size();
//        ref_pos = position;
//
//        // update
//        update();

    return sp2::structure_t();
}
