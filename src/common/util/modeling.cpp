#include <common/neighbor/bond_control_t.hpp>
#include "modeling.hpp"
#include "common/math/vec3_t.hpp"
#include "common/math/vec3_util.hpp"

sp2::structure_t sp2::util::construct_supercell(const sp2::structure_t &input,
    int na, int nb, int nc)
{
    int n_rep = na * nb * nc;

    structure_t supercell = input;

    // lengthen lattice vectors
    int dim[3] = {na, nb, nc};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            supercell.lattice[i][j] *= dim[i];

    // get repeated types, starts at 1 because types is not empty
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

    const auto original_pos = sp2::dtov3(input.positions);
    auto new_pos = std::vector<vec3_t>();

    new_pos.reserve(n_rep * original_pos.size());

    auto add_new = [&](int i, int j, int k) {
        for (std::size_t m = 0; m < original_pos.size(); ++m)
            new_pos.push_back(original_pos[m]
                + i * lattice[0]
                + j * lattice[1]
                + k * lattice[2]);
    };

    // starts at 0 because new_pos is empty
    for (int i = 0; i < na; ++i)
        for (int j = 0; j < nb; ++j)
            for (int k = 0; k < nc; ++k)
                add_new(i, j, k);


    supercell.positions = sp2::v3tod(new_pos);

    return supercell;
}

sp2::structure_t sp2::util::graphene_unit_cell()
{
    // output structure
    structure_t graphene;

    // lattice vectors
    constexpr double lattice_spacing = 2.456047;

    const vec3_t lattice[3] = {
        lattice_spacing * vec3_t{1, 0, 0}, // a1
        lattice_spacing * vec3_t{0.5, std::sqrt(3) / 2, 0}, // a2
        {0, 0, 15} // c, chosen so no vdW interactions
    };

    // copy the lattice
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            graphene.lattice[i][j] = lattice[i][j];

    // atoms
    constexpr double bond_distance = 1.418;
    const vec3_t pos[2] = {
        vec3_t{0, 0, 0},
        vec3_t{std::sqrt(3) / 2, 0.5, 0} * bond_distance
    };

    // set the types
    graphene.types.resize(2, atom_type::CARBON);

    // insert the positions
    for (auto atom : pos)
        for (auto coord : atom)
            graphene.positions.push_back(coord);

    return graphene;
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

sp2::structure_t construct_gnr_impl(bool zigzag, int width, int length,
    bool periodic)
{
    using namespace sp2;

    constexpr double vacuum_sep = 15;
    double lattice[3][3] = {
        {4.254, 0, 0},
        {0, 2.456, 0},
        {0, 0, 0}
    };

    if (zigzag)
        std::swap(lattice[0][0], lattice[1][1]);

    std::vector<atom_type> types(4, atom_type::CARBON);
    std::vector<vec3_t> pos;

    if (zigzag)
    {
        pos = {
            vec3_t{1.228024, 0.000, 0},
            vec3_t{0.000000, 0.709, 0},
            vec3_t{0.000000, 2.127, 0},
            vec3_t{1.228024, 2.836, 0}
        };
    }
    else
    {
        pos = {
            vec3_t{0.709000, 0.000000, 0},
            vec3_t{2.127000, 0.000000, 0},
            vec3_t{0.000000, 1.228024, 0},
            vec3_t{2.836000, 1.228024, 0}
        };
    }

    structure_t unit_cell;
    std::copy_n(lattice[0], 9, unit_cell.lattice[0]);
    unit_cell.types = types;
    unit_cell.positions = sp2::v3tod(pos);

    // repeat the cell in the x direction (width)
    auto supercell = sp2::util::construct_supercell(unit_cell,
        1, (width + 1) / 2);

    // remove extra atoms depending on width
    if (width % 2 == 1)
    {
        supercell.types.resize(supercell.types.size() - 2);
        supercell.positions.resize(supercell.positions.size() - 6);

        // update the lattice (note: technically would make atoms overlap in
        // the armchair gnr case, but since we have vacuum separation it
        // doesn't matter)
        supercell.lattice[1][1] -= lattice[1][1] / 2;
    }

    // repeat the cell in the x direction (length)
    supercell = sp2::util::construct_supercell(supercell, length);

    // add vacuum separation to (x/)y/z
    for (int i = periodic ? 1 : 0; i < 3; ++i)
        supercell.lattice[i][i] += vacuum_sep;

    // center the gnr in the cell
    auto v3pos = sp2::dtov3(supercell.positions);
    vec3_t avg_pos = {};
    for (auto &v : v3pos)
        avg_pos += v;

    avg_pos /= v3pos.size();

    // (max - min) / 2
    auto &sl = supercell.lattice;
    auto offset = vec3_t(sl[0][0], sl[1][1], sl[2][2]) / 2 - avg_pos;

    for (auto &v : v3pos)
        v += offset;

    supercell.positions = sp2::v3tod(v3pos);
    return supercell;
}

sp2::structure_t sp2::util::construct_zz_gnr(int width, int length,
    bool periodic)
{
    return construct_gnr_impl(true, width, length, periodic);
}

sp2::structure_t sp2::util::construct_arm_gnr(int width, int length,
    bool periodic)
{
    return construct_gnr_impl(false, width, length, periodic);
}
