#include <common/neighbor/bond_control_t.hpp>
#include "modeling.hpp"
#include "common/math/vec3_t.hpp"
#include "common/math/vec3_util.hpp"

#include <airebo/system_control_t.hpp>
#include <lammps/lammps_interface.hpp>
#include <common/minimize/minimize.hpp>
#include <common/io/structure.hpp>
#include <run/run_settings_t.hpp>
#include <common/math/rotations.hpp>

void make_agnr_set()
{
    using namespace sp2;

    auto generate = [](int width, int length) {
        std::string filename = std::to_string(width) + "agnr";
        bool periodic = (length == 1);
        if (!periodic)
            filename += "_l" + std::to_string(length);

        filename += ".vasp";

        auto mset = minimize::acgsd_settings_t{};
        mset.output_level = 0;

        sp2::structure_t supercell;
        {
            airebo::system_control_t sys;
            sys.init(util::construct_arm_gnr(width, length, periodic));
            sys.add_hydrogen();

            sys.set_position(
                minimize::acgsd(sys.get_diff_fn(), sys.get_position(), mset)
            );
            supercell = util::construct_supercell(sys.get_structure(), 5, 1, 1);
        }

        auto lset = lammps::lammps_settings_t{};
        lset.compute_lj = false;
        lammps::system_control_t sys(supercell, lset);

        sys.set_position(
            minimize::acgsd(sys.get_diff_fn(), sys.get_position(), mset)
        );
        auto structure = sys.get_structure();
        auto pos = structure.positions;

        double min = std::numeric_limits<double>::max(), max = std::numeric_limits<double>::lowest();
        for (std::size_t i = 0; i < structure.types.size(); ++i)
        {
            if (structure.types[i] == atom_type::HYDROGEN)
                continue;

            min = std::min(min, pos[i].y());
            max = std::max(max, pos[i].y());
        }

        std::cout << width << " " << (max - min) << std::endl;
    };

    for (int width = 3; width < 128; ++width)
    {
        generate(width, 1);

        for (int length = width / 2; length < width * 5; ++length)
            generate(width, length);
    }
}

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

    const auto original_pos = input.positions;
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


    supercell.positions = new_pos;

    return supercell;
}

sp2::structure_t sp2::util::deconstruct_supercell(
    const sp2::structure_t &supercell, int na, int nb, int nc)
{
    int n_rep = na * nb * nc;

    structure_t output = supercell;

    // shorten lattice vectors
    int dim[3] = {na, nb, nc};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            output.lattice[i][j] /= dim[i];

    // number of (actual) atoms is just the size of
    // the supercell / number of cells
    output.types.resize(supercell.types.size() / n_rep);

    // average positions
    vec3_t lattice[3] = {
        vec3_t(output.lattice[0]),
        vec3_t(output.lattice[1]),
        vec3_t(output.lattice[2])
    };

    auto &pos = output.positions;
    pos.resize(output.types.size());
    std::fill(pos.begin(), pos.end(), sp2::vec3_t{0, 0, 0});

    int counter = 0;
    auto avg_old = [&](int i, int j, int k) {
        for (std::size_t m = 0; m < pos.size(); ++m)
            pos[m] += (supercell.positions[counter * pos.size() + m]
                - i * lattice[0]
                - j * lattice[1]
                - k * lattice[2]) / n_rep;

        counter++;
    };

    // starts at 0 because new_pos is empty
    for (int i = 0; i < na; ++i)
        for (int j = 0; j < nb; ++j)
            for (int k = 0; k < nc; ++k)
                avg_old(i, j, k);


    return output;
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

    // atom types
    graphene.types = {atom_type::CARBON, atom_type::CARBON};

    // atom positions
    constexpr double bond_distance = 1.418;
    graphene.positions = {
        vec3_t{0, 0, 0},
        vec3_t{std::sqrt(3) / 2, 0.5, 0} * bond_distance
    };

    return graphene;
}

sp2::structure_t sp2::util::make_hydrogen_terminated(
    const sp2::structure_t &input)
{
    // maximum distance between two atoms to be considered bonded in Angstroms
    constexpr double max_dist = 1.6;

    // use fbc::bond_control_t to calculate all the bonds between atoms
    // in the system
    fbc::bond_control_t bond_control;
    bond_control.init(input.lattice, max_dist, 0);
    bond_control.update(sp2::v3tod(input.positions));

    auto types = input.types;
    auto pos = input.positions;

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

//        vec3_t bonds[3];
//        int neigh_ids[3];
//
//        for (int j = 0; j < n_bonds; ++j)
//        {
//            auto edge = graph.edge(idx, j);
//            bonds[j] = bond_deltas[edge.id];
//            neigh_ids[j] = edge.b;
//        }
//
//        // initial CH bond length for added atoms
//        constexpr double ch_dist = 1.09;
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
    unit_cell.positions = pos;

    // repeat the cell in the x direction (width)
    auto supercell = sp2::util::construct_supercell(unit_cell,
        1, (width + 1) / 2);

    // remove extra atoms depending on width
    if (width % 2 == 1)
    {
        supercell.types.resize(supercell.types.size() - 2);
        supercell.positions.resize(supercell.positions.size() - 2);

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
    auto v3pos = supercell.positions;
    vec3_t avg_pos = {};
    for (auto &v : v3pos)
        avg_pos += v;

    avg_pos /= v3pos.size();

    // (max - min) / 2
    auto &sl = supercell.lattice;
    auto offset = vec3_t(sl[0][0], sl[1][1], sl[2][2]) / 2 - avg_pos;

    for (auto &v : v3pos)
        v += offset;

    supercell.positions = v3pos;
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

sp2::structure_t sp2::util::center_by_avg(const sp2::structure_t &input)
{
    auto pos = input.positions;

    vec3_t avg_pos = {};
    for (auto &atom : pos)
        avg_pos += atom;

    avg_pos /= pos.size();

    vec3_t lattice_center = (
        vec3_t(input.lattice[0]) +
        vec3_t(input.lattice[1]) +
        vec3_t(input.lattice[2])
    ) / 2;

    auto delta = lattice_center - avg_pos;
    for (auto &atom : pos)
        atom += delta;

    auto new_structure = input;
    new_structure.positions = pos;

    return new_structure;
}

sp2::structure_t sp2::util::make_capped_nanotube(const sp2::structure_t &cap,
    const sp2::structure_t &segment, int n_segment)
{
    structure_t output = cap,
        rep = util::construct_supercell(segment, 1, 1, n_segment);

    output.positions.insert(output.positions.end(),
        rep.positions.begin(), rep.positions.end());

    output.types.insert(output.types.end(),
        rep.types.begin(), rep.types.end());

    return output;
}


sp2::vec3_t get_tvec(const sp2::structure_t &input)
{
    auto temp = input.positions;
    std::sort(temp.begin(), temp.end(), [](auto &a, auto &b) {
        return a.z() < b.z();
    });

    sp2::vec3_t top_avg = {};
    for (int i = 0; i < 6; ++i)
        top_avg += temp[temp.size() - i - 1];

    sp2::vec3_t bottom_avg = {};
    for (int i = 0; i < 12; ++i)
        bottom_avg += temp[i];

    return (top_avg / 6 - bottom_avg / 12).unit_vector();
}

void shorten_tube(sp2::structure_t &input, double delta)
{
    using namespace sp2;

    auto &pos = input.positions;

    // sort by z
    std::sort(pos.begin(), pos.end(),
        [](auto &a, auto &b) {
            return a.z() < b.z();
    });

    auto get_delta = [bottom = pos.back().z(), &pos]() -> double {
        return pos.back().z() - bottom;
    };

    double lattice[3][3] = {};
    sp2::fbc::bond_control_t bc;
    bc.init(lattice, 1.7, 0);

    std::vector<int> n_bonds;
    auto update_bonds = [&]{
        bc.update(sp2::v3tod(pos));
        auto graph = bc.get_graph();

        n_bonds.clear();
        for (std::size_t i = 0; i < graph.n_vertices(); ++i)
            n_bonds.push_back(graph.degree(i));
    };

    while (!pos.empty() && get_delta() < delta)
    {
        pos.pop_back();
        update_bonds();

        for (std::size_t i = 0; i < pos.size();)
        {
            if (n_bonds[i] > 1)
                i += 1;
            else
            {
                pos.erase(pos.begin() + i);
                n_bonds.erase(n_bonds.begin() + i);
                i -= 1;
            }
        }
    }

    input.types.resize(input.positions.size());
}

void make_dataset(const sp2::run_settings_t &, MPI_Comm)
{
    using namespace sp2;

    structure_t tube;
    if (!io::read_structure("tube5-5.xyz", tube))
        throw std::runtime_error("cant open file");


}
