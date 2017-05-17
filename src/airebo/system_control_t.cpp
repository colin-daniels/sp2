#include "airebo/system_control_t.hpp"
#include "airebo/utility_functions.hpp"
#include "airebo/airebo_util.hpp"
#include "airebo/interpolation_coeff.hpp"

#include "common/io/util.hpp"
#include "common/io/structure.hpp"
#include "common/math/rotations.hpp"
#include "common/structure_t.hpp"
#include "common/math/vec3_t.hpp"
#include "common/math/vec3_util.hpp"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <queue>

using namespace std;
using namespace sp2;
using namespace sp2::airebo;

// initialization function, takes lattice vectors, initial atom positions, and types
void airebo::system_control_t::init(const double lattice[3][3],
    const vector<double> &position_in,
    const vector<atom_type> &types_in)
{
    copy_n(lattice[0], 9, ref_lattice[0]);
    bond_control.init(lattice,
        D_max[static_cast<int>(bond_type::CARBON_CARBON)], 0.1);
    position = position_in;
    types = types_in;
    force = vector<double>(position.size(), 0);
    static_atom = vector<bool>(position.size() / 3, false);
    update();
}

void airebo::system_control_t::init(const structure_t &structure)
{
    init(structure.lattice, sp2::v3tod(structure.positions), structure.types);
}

structure_t airebo::system_control_t::get_structure() const
{
    structure_t structure;

    structure.types = types;
    structure.positions = sp2::dtov3(get_position());
    copy_n(ref_lattice[0], 9, structure.lattice[0]);

    return structure;
}

// wrap positions into the current unit cell
void airebo::system_control_t::wrap_positions()
{
    bond_control.wrap_positions(position);
}

// set/update lattice, note: automatically wraps positions into unit cell before modifying the lattice
void airebo::system_control_t::set_lattice(const double lattice[3][3])
{
    for (int i = 0; i < 9; ++i)
    {
        if (lattice[i/3][i%3] != ref_lattice[i/3][i%3])
        {
            copy_n(lattice[0], 9, ref_lattice[0]);
            bond_control.set_lattice(lattice);
            break;
        }
    }
}

void airebo::system_control_t::set_structure(const structure_t &input)
{
    types = input.types;
    position = sp2::v3tod(input.positions);
    for (int i = 0; i < 9; ++i)
    {
        if (input.lattice[i/3][i%3] != ref_lattice[i/3][i%3])
        {
            copy_n(input.lattice[0], 9, ref_lattice[0]);
            bond_control.set_lattice(input.lattice);
            break;
        }
    }
}

void airebo::system_control_t::write_output(std::string filename)
{
    io::write_structure(filename, get_structure());
}

void airebo::system_control_t::append_output(std::string filename)
{
    io::write_structure(filename, get_structure(), true);
}

// update function, calculates total potential and forces for the current atom positions
void airebo::system_control_t::update()
{
    na = position.size() / 3;
    static_atom.resize(na, false);

    potential.resize(na, 0);
    fill(potential.begin(), potential.end(), 0.0);

    if ((int)types.size() != na)
        types.resize(na, atom_type::CARBON);

    if (ref_pos.size() > position.size())
        ref_pos.resize(position.size());
    else if (ref_pos.size() < position.size())
        for (size_t i = ref_pos.size(); i < position.size(); ++i)
            ref_pos.push_back(position[i]);

    for (size_t i = 0; i < position.size(); ++i)
        if (static_atom[i / 3])
            position[i] = ref_pos[i];

    // update bonds
    bond_control.update(position);

    // update number of bonds
    nb = bond_control.get_nb();

    // reset variables
    total_potential = 0;

    force.resize(na * 3);
    fill(force.begin(), force.end(), 0);

    // stage 1
    update_stage1();

    // stage 2
    update_stage2();

    // stage 3
    update_stage3();

    // stage 4
    update_stage4();

    // stage 5
    update_stage5();

    // stage 6
    update_stage6();

    // stage 7
    update_stage7();

    // distribute forces to atoms
    update_distribute();
}

void airebo::system_control_t::update_stage1()
{
    // reset variables
    cutoff.resize(nb * 2);
    attr_p.resize(nb * 2);
    bond_force.resize(nb * 3);
    bond_dir_force.resize(nb);
    Na_vars.resize(na * 4);
    b_types.resize(nb);

    fill(cutoff.begin(), cutoff.end(), 0);
    fill(attr_p.begin(), attr_p.end(), 0);
    fill(bond_force.begin(), bond_force.end(), 0);
    fill(bond_dir_force.begin(), bond_dir_force.end(), 0);
    fill(Na_vars.begin(), Na_vars.end(), 0);
    fill(b_types.begin(), b_types.end(), bond_type::CARBON_CARBON);

    // references
    const vector<double> &lengths = bond_control.get_bond_lengths();
    const vector<int> &bond_ids = bond_control.get_bond_ids(),
        &offsets = bond_control.get_bond_offsets(),
        &sister_ids = bond_control.get_sister_bonds();

    for (int i = 0; i < na; ++i)
    {
        for (int j = offsets[i]; j < offsets[i + 1]; ++j)
        {
            // build the bond type list
            int id_b = bond_ids[j];
            b_types[j] = get_bond_type(types[i], types[id_b]);

            if (i < id_b)
                continue;

            // get the cutoff function
            cutoff_fn(lengths[j], (int)b_types[j], cutoff[j * 2],
                cutoff[j * 2 + 1]);

            if (cutoff[j * 2] == 0)
                continue;

            // start calculating the N_variables, get N_t and N_h/N_c
            Na_vars[   i * 4] += cutoff[j * 2];
            Na_vars[id_b * 4] += cutoff[j * 2];
            Na_vars[i * 4 + 1 + (int)types[id_b]] += cutoff[j * 2];
            Na_vars[id_b * 4 + 1 + (int)types[i]] += cutoff[j * 2];

            // get the attractive and repulsive potentials
            double rep_p,
                d_rep_p;
            lk_ar_ptnl(lengths[j], (int)b_types[j], rep_p, d_rep_p,
                attr_p[j * 2], attr_p[j * 2 + 1]);

            // add the repulsive potential to the total potential now, since no other terms are needed
            total_potential += cutoff[j * 2] * rep_p;
            bond_dir_force[j] += cutoff[j * 2 + 1] * rep_p + cutoff[j * 2] * d_rep_p;

            potential[i] += cutoff[j * 2] * rep_p;
            potential[id_b] += cutoff[j * 2] * rep_p;


            // set the values for the sister bond (id_a < id_b)
            for (int k = 0; k < 2; ++k)
            {
                cutoff[sister_ids[j] * 2 + k] = cutoff[j * 2 + k];
                attr_p[sister_ids[j] * 2 + k] = attr_p[j * 2 + k];
            }
        }
    }

    // Finish getting Na_vars (half of Na_conj)
    for (int i = 0; i < na; ++i)
        for (int j = offsets[i]; j < offsets[i + 1]; ++j)
            if (cutoff[j * 2] != 0 && types[bond_ids[j]] == atom_type::CARBON)
                Na_vars[i * 4 + 3] += cutoff[j * 2] * F_conj(Na_vars[bond_ids[j] * 4] - cutoff[j * 2]);
}

void airebo::system_control_t::update_stage2()
{
    // reset variables
    Nb_vars.resize(nb * 7);
    fill(Nb_vars.begin(), Nb_vars.end(), 0);

    // references
    const vector<int> &bond_ids = bond_control.get_bond_ids(),
        &offsets = bond_control.get_bond_offsets(),
        &sister_ids = bond_control.get_sister_bonds();

    // calculate N variables for bonds
    for (int i = 0; i < na; ++i)
    {
        for (int j = offsets[i]; j < offsets[i + 1]; ++j)
        {
            if (cutoff[j * 2] == 0)
                continue;

            // N_ti for this bond and N_tj for the sister bond
            Nb_vars[j * 7] = Na_vars[i * 4] - cutoff[j * 2];
            Nb_vars[sister_ids[j] * 7 + 1] = Nb_vars[j * 7];

            // N_h, N_c
            Nb_vars[j * 7 + 2] = Na_vars[i * 4 + 1];
            Nb_vars[j * 7 + 3] = Na_vars[i * 4 + 2];
            Nb_vars[j * 7 + 2 + (int)types[bond_ids[j]]] -= cutoff[j * 2];

            // Nca for j, Ncb for sister[j]
            Nb_vars[j * 7 + 5] = Na_vars[i * 4 + 3];
            if (types[bond_ids[j]] == atom_type::CARBON)
            {
                Nb_vars[j * 7 + 5] -= cutoff[j * 2]
                                      * F_conj(Na_vars[bond_ids[j] * 4] - cutoff[j * 2]);
            }
            Nb_vars[sister_ids[j] * 7 + 6] = Nb_vars[j * 7 + 5];
        }
    }

    // do the N_conj sums
    for (int i = 0; i < nb; ++i)
    {
        if (cutoff[i * 2] == 0)
            continue;

        Nb_vars[i * 7 + 4] =  min(1 + Nb_vars[i * 7 + 5] * Nb_vars[i * 7 + 5]
                                  + Nb_vars[i * 7 + 6] * Nb_vars[i * 7 + 6], 9.0);
        Nb_vars[i * 7 + 5] *= 2;
        Nb_vars[i * 7 + 6] *= 2;
    }
}

void airebo::system_control_t::update_vdw()
{
    if (!vdw_enabled)
        return;

    const vector<int> &offsets = bond_control.get_bond_offsets(),
        &bond_ids = bond_control.get_bond_ids(),
        &sister_ids = bond_control.get_sister_bonds();

    vdw_force.resize(na * 3, 0);
    bond_control.calc_delta_all(vdw_deltas, vdw_working);

    vector<int> cc_paths(na, 0);
    vector<double> temp_visited(na, 0);
    for (int i = 0, total = 0; i < na; ++i)
    {
        // connectivity term
        fill(temp_visited.begin(), temp_visited.end(), 0);
        fill(cc_paths.begin(), cc_paths.end(), 0);

        queue<int> to_visit;
        to_visit.push(i);
        temp_visited[i] = 1;

        // breadth first search, keep maximum connection path
        for (int j = 0; j < 3; ++j)
        {
            int current_size = to_visit.size();
            for (int k = 0; k < current_size; ++k)
            {
                int id = to_visit.front();
                for (int l = offsets[id]; l < offsets[id + 1]; ++l)
                {
                    int id_b = bond_ids[l];
                    to_visit.pop();

                    if (temp_visited[id_b] == 0)
                        to_visit.push(id_b);

                    if (temp_visited[id_b] < cutoff[l * 2] * temp_visited[id])
                    {
                        temp_visited[id_b] = cutoff[l * 2] * temp_visited[id];
                        cc_paths[id_b] = sister_ids[l];
                    }
                }
            }
        }

        // get actual forces
        for (int j = i + 1; j < na; ++j, ++total)
        {
            // don't calculate if third nearest neighbor or less
            if (temp_visited[j] == 1)
                continue;

            int btype = (int)get_bond_type(types[i], types[j]);
            double cc_val = (1 - temp_visited[j]);

            double *delta = &vdw_deltas[total * 3];

            double len = sqrt(delta[0] * delta[0] +
                              delta[1] * delta[1] +
                              delta[2] * delta[2]);

            // calculate the Lennard-Jones potential and associated force
            double lj_val, lj_force;
            lj_12_6(len, btype, lj_val, lj_force);

            double S_r = 0, dS_r = 0;
            if (len < lj_D_min[btype])
                S_r = 1;
            else if (len < lj_D_max[btype])
            {
                lj_S_r(len, btype, S_r, dS_r);

                // b_ij = F_ij + T_ij + 0.5 * b_ij^sigma

            }

            // distribute forces on partial connectivity
            if (cc_val < 1)
                for (int k = cc_paths[j]; k != 0; k = cc_paths[bond_ids[k]])
                    bond_dir_force[k] += (1 - S_r) * -temp_visited[j] * lj_val
                                         * cutoff[k * 2 + 1] / cutoff[k * 2];

            // update force and potential
            lj_force = (1 - S_r) * cc_val * lj_force
                       + -dS_r * cc_val * lj_val / len;
            lj_val = (1 - S_r) * cc_val * lj_val;

            // add to the total potential
            total_potential += lj_val;

            // distribute forces in the bond direction
            for (int k = 0; k < 3; ++k)
            {
                force[i * 3 + k] += lj_force * delta[k];
                force[j * 3 + k] -= lj_force * delta[k];
            }
        }
    }
}

void airebo::system_control_t::update_stage3()
{
    // reset variables
    F.resize(nb * 4);
    T.resize(nb * 4);
    P.resize(nb * 3);

    fill(F.begin(), F.end(), 0);
    fill(T.begin(), T.end(), 0);
    fill(P.begin(), P.end(), 0);

    // references
    const vector<int> &bond_ids = bond_control.get_bond_ids(),
        &offsets = bond_control.get_bond_offsets();

    for (int i = 0; i < na; ++i)
    {
        // get P values
        for (int j = offsets[i]; j < offsets[i + 1]; ++j)
        {
            if (cutoff[j * 2] == 0)
                continue;

            // get the interpolation id
            double N_h = Nb_vars[j * 7 + 2],
                N_c = Nb_vars[j * 7 + 3];
            int id_p = ((int)N_h) * 5 + (int)N_c;

            if (N_h >= 4 || N_c >= 4)
                continue;

            if (N_c == std::floor(N_c) && N_h == std::floor(N_h))
            {
                if (b_types[j] == bond_type::CARBON_CARBON)
                    P[j * 3] = PCC[(int)N_h][(int)N_c];
                else if (b_types[j] == bond_type::HYDROGEN_CARBON && types[i] == atom_type::CARBON)
                    P[j * 3] = PCH[(int)N_h][(int)N_c];
                continue;
            }

            // interpolation coefficient pointer
            const double *coeff_ptr = nullptr;

            if (b_types[j] == bond_type::CARBON_CARBON)
                coeff_ptr = map_pcc[id_p] != -1 ? coeff_pcc[map_pcc[id_p]] : nullptr;
            else if (b_types[j] == bond_type::HYDROGEN_CARBON && types[i] == atom_type::CARBON)
                coeff_ptr = map_pch[id_p] != -1 ? coeff_pch[map_pch[id_p]] : nullptr;

            // dont calculate if the value of P is zero
            if (coeff_ptr == nullptr)
                continue;

            // get the inputs for the bicubic spline
            double inp[4][4] = {
                {1, N_h, N_h * N_h, N_h * N_h * N_h},
                {0,   1,   2 * N_h,   3 * N_h * N_h},
                {1, N_c, N_c * N_c, N_c * N_c * N_c},
                {0,   1,   2 * N_c,   3 * N_c * N_c}
            };

            // xy, dxy, xdy
            double input[48] = {};
            loop_dger<8, 4>(inp[0], inp[2],  &input[0]);
            loop_dger<4, 4>(inp[0], inp[3], &input[32]);

            // get P[3] = {xy, dxy, xdy}
            loop_dgemv<3, 16>(input, coeff_ptr, &P[j * 3]);
        }

        // get T and F values
        for (int j = offsets[i]; j < offsets[i + 1]; ++j)
        {
            if (i < bond_ids[j] || cutoff[j * 2] == 0)
                continue;

            // get the interpolation id
            double N_ti = min(Nb_vars[j * 7], 3.0),
                N_tj = min(Nb_vars[j * 7 + 1], 3.0),
                N_conj = min(Nb_vars[j * 7 + 4], 9.0);
            int id_T = ((int)N_ti) * 40 + ((int)N_tj) * 10 + (int)N_conj,
                id_F = id_T;

            if (N_ti == std::floor(N_ti) &&
                N_tj == std::floor(N_tj) &&
                N_conj == std::floor(N_conj))
            {
                if (b_types[j] == bond_type::CARBON_CARBON)
                {
                    T[j * 4] = TCC[(int)N_ti][(int)N_tj][(int)N_conj];
                    F[j * 4] = FCC[(int)N_ti][(int)N_tj][(int)N_conj];
                    F[j * 4 + 1] = DXFCC[(int)N_ti][(int)N_tj][(int)N_conj];
                    F[j * 4 + 2] = DYFCC[(int)N_ti][(int)N_tj][(int)N_conj];
                    F[j * 4 + 3] = DZFCC[(int)N_ti][(int)N_tj][(int)N_conj];
                }
                else if (b_types[j] == bond_type::HYDROGEN_CARBON)
                    F[j * 4] = FCH[(int)N_ti][(int)N_tj][(int)N_conj];
                else
                    F[j * 4] = FHH[(int)N_ti][(int)N_tj][(int)N_conj];
                continue;
            }

            // coefficient pointers
            const double *t_coeff = nullptr,
                *f_coeff = nullptr;

            if (b_types[j] == bond_type::CARBON_CARBON)
            {
                t_coeff = map_tcc[id_T] != -1 ? coeff_tcc[map_tcc[id_T]] : nullptr;
                f_coeff = map_fcc[id_F] != -1 ? coeff_fcc[map_fcc[id_F]] : nullptr;
            }
            else if (b_types[j] == bond_type::HYDROGEN_CARBON)
                f_coeff = map_fch[id_F] != -1 ? coeff_fch[map_fch[id_F]] : nullptr;
            else
                f_coeff = map_fhh[id_F] != -1 ? coeff_fhh[map_fhh[id_F]] : nullptr;

            // only keep calculating if one of the quantities needs to be evaluated
            if (t_coeff == nullptr && f_coeff == nullptr)
                continue;

            // get the inputs for the tricubic spline
            double inp[6][4] = {
                {1,   N_ti,     N_ti * N_ti,       N_ti * N_ti * N_ti},
                {0,      1,        2 * N_ti,          3 * N_ti * N_ti},
                {1,   N_tj,     N_tj * N_tj,       N_tj * N_tj * N_tj},
                {0,      1,        2 * N_tj,          3 * N_tj * N_tj},
                {1, N_conj, N_conj * N_conj, N_conj * N_conj * N_conj},
                {0,      1,      2 * N_conj,      3 * N_conj * N_conj}
            };

            // yz dyz ydz
            double working[48] = {};
            loop_dger<8, 4>(inp[2], inp[4],  &working[0]);
            loop_dger<4, 4>(inp[2], inp[5], &working[32]);

            // xyz dxyz xdyz xydz
            double input[256] = {};
            loop_dger<4, 16>(inp[0],  &working[0],   &input[0]);
            loop_dger<4, 16>(inp[1],  &working[0],  &input[64]);
            loop_dger<4, 16>(inp[0], &working[16], &input[128]);
            loop_dger<4, 16>(inp[0], &working[32], &input[192]);

            // get T[4] = {xyz, dxyz, xdyz, xydz}
            if (t_coeff != nullptr)
                loop_dgemv<4, 64>(input, t_coeff, &T[j * 4]);

            // get F[4] = {xyz, dxyz, xdyz, xydz}
            if (f_coeff != nullptr)
                loop_dgemv<4, 64>(input, f_coeff, &F[j * 4]);
        }
    }
}

void airebo::system_control_t::update_stage4()
{
    // references
    const vector<double> &lengths = bond_control.get_bond_lengths(),
        &deltas = bond_control.get_bond_deltas();
    const vector<int> &bond_ids = bond_control.get_bond_ids(),
        &offsets = bond_control.get_bond_offsets();

    // get offsets
    sd_offsets.resize(na + 1);
    od_offsets.resize(nb + 1);

    sd_offsets[0] = 0;
    od_offsets[0] = 0;

    for (int i = 0; i < na; ++i)
    {
        // number of bond for the atom i
        int nb_a = offsets[i + 1] - offsets[i];

        // 'self' dot product offset
        sd_offsets[i + 1] = sd_offsets[i] + nb_a * nb_a;

        // 'other'/'inter atom' dot product offset
        for (int j = offsets[i]; j < offsets[i + 1]; ++j)
        {
            int nb_b = offsets[bond_ids[j] + 1] - offsets[bond_ids[j]];

            // these don't need to be calculated if the following conditions are true
            if (i < bond_ids[j] || T[j * 4] == 0 || cutoff[j * 2] == 0)
                od_offsets[j + 1] = od_offsets[j];
            else
                od_offsets[j + 1] = od_offsets[j] + nb_b * nb_a;
        }
    }

    // allocate vectors now that the needed space has been calculated
    other_dots.resize(od_offsets.back());
    self_dots.resize(sd_offsets.back());
    cos_vals.resize(sd_offsets.back());

    // get the 'self' dot products, for G(theta), and the 'inter-atom'/'other' ones for the dihedral force term
    for (int i = 0; i < na; ++i)
    {
        // number of bond for the atom i
        const int nb_a = offsets[i + 1] - offsets[i];
        if (nb_a == 0)
            continue;

        const double *delta_ptr_a = &deltas[offsets[i] * 3],
            *len_ptr_a   = &lengths[offsets[i]];

        double *sd_ptr  = &self_dots[sd_offsets[i]],
            *cos_ptr = &cos_vals[sd_offsets[i]];

        // calculate dot products and cosines between bonds
        for (int j = 0; j < nb_a; ++j)
        {
            sd_ptr[j * nb_a + j] = len_ptr_a[j] * len_ptr_a[j];

            for (int k = j + 1; k < nb_a; ++k)
            {
                double val = delta_ptr_a[j * 3] * delta_ptr_a[k * 3] +
                             delta_ptr_a[j * 3 + 1] * delta_ptr_a[k * 3 + 1] +
                             delta_ptr_a[j * 3 + 2] * delta_ptr_a[k * 3 + 2];

                // normal dot products
                sd_ptr[j * nb_a + k]  = val;
                sd_ptr[k * nb_a + j]  = val;

                // cosine values
                cos_ptr[j * nb_a + k] = val / (len_ptr_a[j] * len_ptr_a[k]);
                cos_ptr[k * nb_a + j] = val / (len_ptr_a[j] * len_ptr_a[k]);

            }
        }

        // get the inter atom dot products
        for (int j = offsets[i]; j < offsets[i + 1]; ++j)
        {
            // continue if no space was allocated for this calculation (aka, it isn't needed)
            if (od_offsets[j + 1] == od_offsets[j])
                continue;

            // number of bonds for the other atom
            const int nb_b = offsets[bond_ids[j] + 1] - offsets[bond_ids[j]];

            double *od_ptr = &other_dots[od_offsets[j]];
            const double *delta_ptr_b = &deltas[offsets[bond_ids[j]] * 3];

            for (int k = 0; k < nb_a; ++k)
            {
                for (int p = 0; p < nb_b; ++p)
                {
                    od_ptr[k * nb_b + p] = delta_ptr_a[k * 3] * delta_ptr_b[p * 3] +
                                           delta_ptr_a[k * 3 + 1] * delta_ptr_b[p * 3 + 1] +
                                           delta_ptr_a[k * 3 + 2] * delta_ptr_b[p * 3 + 2];
                }
            }
        }
    }
}

void airebo::system_control_t::update_stage5()
{
    // reset variables
    force_on_cutoff.resize(na);
    fill(force_on_cutoff.begin(), force_on_cutoff.end(), 0);

    // references
    const vector<double> &lengths = bond_control.get_bond_lengths(),
        &deltas = bond_control.get_bond_deltas();
    const vector<int> &bond_ids = bond_control.get_bond_ids(),
        &offsets = bond_control.get_bond_offsets();

    // get the angular part (gtheta sum) of the bond order term
    for (int i = 0; i < na; ++i)
    {
        int nb_a = offsets[i + 1] - offsets[i];
        if (nb_a == 0)
            continue;

        double *cos_ptr = &cos_vals[sd_offsets[i]];

        vector<double> gtheta(nb_a),
            dgtheta(nb_a * 2);

        for (int j = offsets[i], j2 = 0; j < offsets[i + 1]; ++j, ++j2)
        {
            if (cutoff[j * 2] == 0)
                continue;

            fill(gtheta.begin(), gtheta.end(), 0);
            fill(dgtheta.begin(), dgtheta.end(), 0);

            double b_sum = 1 + P[j * 3];
            for (int k = offsets[i], k2 = 0; k < offsets[i + 1]; ++k, ++k2)
            {
                if (cutoff[k * 2] == 0 || k == j)
                    continue;

                if (types[i] == atom_type::CARBON)
                    gtheta_C(cos_ptr[j2 * nb_a + k2], Nb_vars[j * 7],
                        gtheta[k2], dgtheta[k2 * 2], dgtheta[k2 * 2 + 1]);
                else
                {
                    gtheta_H(cos_ptr[j2 * nb_a + k2], gtheta[k2], dgtheta[k2 * 2]);
                    if (b_types[j] == bond_type::HYDROGEN_HYDROGEN &&
                        b_types[k] == bond_type::HYDROGEN_HYDROGEN)
                        gtheta[k2] *= exp(4.0);
                }

                b_sum += gtheta[k2] * cutoff[k * 2];
            }

            b_sum = 1.0 / sqrt(b_sum);

            double val = cutoff[j * 2] * -0.5 * attr_p[j * 2] * b_sum;
            total_potential += val;

            potential[i] += val;
            potential[bond_ids[j]] += val;

            bond_dir_force[j] += -0.5 * b_sum * (cutoff[j * 2 + 1] * attr_p[j * 2]
                                                 + cutoff[j * 2] * attr_p[j * 2 + 1]);
            double mult = 0.25 * cutoff[j * 2] * attr_p[j * 2] * b_sum * b_sum * b_sum;


            for (int k = offsets[i], k2 = 0; k < offsets[i + 1]; ++k, ++k2)
            {
                if (k == j || cutoff[k * 2] == 0)
                    continue;

                bond_dir_force[k] += mult * (cutoff[k * 2] * dgtheta[k2 * 2]
                                             * -cos_ptr[j2 * nb_a + k2] / lengths[k] + cutoff[k * 2 + 1]
                                                                                       * (gtheta[k2] + P[j * 3 + 1 + (int)types[bond_ids[k]]]));
                bond_dir_force[j] += mult * (cutoff[k * 2] * dgtheta[k2 * 2]
                                             * -cos_ptr[j2 * nb_a + k2] / lengths[j]);

                for (int p = 0; p < 3; ++p)
                {
                    bond_force[k * 3 + p] += mult * cutoff[k * 2]
                                             * dgtheta[k2 * 2] * deltas[j * 3 + p] / (lengths[j] * lengths[k]);
                    bond_force[j * 3 + p] += mult * cutoff[k * 2]
                                             * dgtheta[k2 * 2] * deltas[k * 3 + p] / (lengths[j] * lengths[k]);
                }

                if (dgtheta[k2 * 2 + 1] != 0)
                {
                    force_on_cutoff[i] += mult * cutoff[k * 2]
                                          * dgtheta[k2 * 2 + 1];
                    bond_dir_force[j] -= mult * cutoff[k * 2]
                                         * dgtheta[k2 * 2 + 1] * cutoff[j * 2 + 1];
                }
            }
        }
    }
}

void airebo::system_control_t::update_stage6()
{
    // references
    const vector<double> &deltas = bond_control.get_bond_deltas();
    const vector<int> &bond_ids = bond_control.get_bond_ids(),
        &offsets = bond_control.get_bond_offsets(),
        &sister_ids = bond_control.get_sister_bonds();

    vector<double> tvec = self_dots;
    for (int i = 0; i < na; ++i)
    {
        int nb_a = offsets[i + 1] - offsets[i];
        double *tv_ptr = &tvec[sd_offsets[i]]; // array -> [nb_a][nb_a]
        for (int j = 0; j < nb_a; ++j)
            if (types[i] == atom_type::CARBON)
                for (int k = j + 1; k < nb_a; ++k)
                    tv_ptr[k * nb_a + j] = tv_ptr[nb_a * j + j]
                                           * tv_ptr[nb_a * k + k] - tv_ptr[k * nb_a + j] * tv_ptr[k * nb_a + j];
    }

    for (int id_a = 0; id_a < na; ++id_a)
    {
        int nb_a = offsets[id_a + 1] - offsets[id_a];
        double *tv_a = &tvec[sd_offsets[id_a]]; // array -> [nb_a][nb_a]

        for (int i = 0; i < nb_a; ++i)
        {
            int id_b = bond_ids[i + offsets[id_a]],
                nb_b = offsets[id_b + 1] - offsets[id_b];

            if (nb_b == 0 || T[(i + offsets[id_a]) * 4] == 0)
                continue;

            int j = sister_ids[i + offsets[id_a]] - offsets[id_b];
            double *tv_b = &tvec[sd_offsets[id_b]],
                *od_ptr = &other_dots[od_offsets[i + offsets[id_a]]];

            // dihedral sum
            double dih_sum = 0;

            for (int k = 0; k < nb_a; ++k)
            {
                double c_k = cutoff[(k + offsets[id_a]) * 2],
                    dc_k = cutoff[(k + offsets[id_a]) * 2 + 1];
                if (i == k || c_k == 0)
                    continue;

                for (int l = 0; l < nb_b; ++l)
                {
                    double c_l = cutoff[(l + offsets[id_b]) * 2],
                        dc_l = cutoff[(l + offsets[id_b]) * 2 + 1];
                    if (j == l || c_l == 0)
                        continue;

                    double denom_a = tv_a[i * nb_a + k],
                        denom_b = tv_b[j * nb_b + l],
                        dot_a = tv_a[k * nb_a + i],
                        dot_b = tv_b[l * nb_b + j];

                    if (i < k) {swap(denom_a, dot_a);};
                    if (j < l) {swap(denom_b, dot_b);};

                    double numerator = (od_ptr[i * nb_b + j]
                                        * od_ptr[k * nb_b + l] - dot_a * dot_b);

                    double val = (numerator * numerator) / (denom_a * denom_b);
                    dih_sum += c_l * c_k * (1 - val);

                    double mult = cutoff[(i + offsets[id_a]) * 2]
                                  * -attr_p[(i + offsets[id_a]) * 2] * T[(i + offsets[id_a]) * 4];

                    bond_dir_force[k + offsets[id_a]] += mult * c_l * dc_k * (1 - val);
                    bond_dir_force[l + offsets[id_b]] += mult * dc_l * c_k * (1 - val);

                    mult *= c_l * c_k;

                    for (int p = 0; p < 3; ++p)
                    {
                        bond_force[(i + offsets[id_a]) * 3 + p] += mult * val *
                                                                   2 * ((tv_a[k * nb_a + k] * deltas[(i + offsets[id_a]) * 3 + p]
                                                                         - dot_a * deltas[(k + offsets[id_a]) * 3 + p]) / denom_a
                                                                        - (deltas[(j + offsets[id_b]) * 3 + p] * od_ptr[k * nb_b + l]
                                                                           - deltas[(k + offsets[id_a]) * 3 + p] * dot_b) / numerator);
                        bond_force[(k + offsets[id_a]) * 3 + p] += mult * val *
                                                                   2 * ((tv_a[i * nb_a + i] * deltas[(k + offsets[id_a]) * 3 + p]
                                                                         - dot_a * deltas[(i + offsets[id_a]) * 3 + p]) / denom_a
                                                                        - (od_ptr[i * nb_b + j] * deltas[(l + offsets[id_b]) * 3 + p]
                                                                           - deltas[(i + offsets[id_a]) * 3 + p] * dot_b) / numerator);

                        bond_force[(j + offsets[id_b]) * 3 + p] += mult * val *
                                                                   2 * ((tv_b[l * nb_b + l] * deltas[(j + offsets[id_b]) * 3 + p]
                                                                         - dot_b * deltas[(l + offsets[id_b]) * 3 + p]) / denom_b
                                                                        - (deltas[(i + offsets[id_a]) * 3 + p] * od_ptr[k * nb_b + l]
                                                                           - deltas[(l + offsets[id_b]) * 3 + p] * dot_a) / numerator);
                        bond_force[(l + offsets[id_b]) * 3 + p] += mult * val *
                                                                   2 * ((tv_b[j * nb_b + j] * deltas[(l + offsets[id_b]) * 3 + p]
                                                                         - dot_b * deltas[(j + offsets[id_b]) * 3 + p]) / denom_b
                                                                        - (od_ptr[i * nb_b + j] * deltas[(k + offsets[id_a]) * 3 + p]
                                                                           - deltas[(j + offsets[id_b]) * 3 + p] * dot_a) / numerator);
                    }
                }
            }

            for (int k = 0; k < 4; ++k)
            {
                T[(offsets[id_a] + i) * 4 + k] *= dih_sum;
            }
        }
    }
}

void airebo::system_control_t::update_stage7()
{
    // references
    const vector<int> &bond_ids = bond_control.get_bond_ids(),
        &offsets = bond_control.get_bond_offsets(),
        &sister_ids = bond_control.get_sister_bonds();

    for (int i = 0; i < na; ++i)
    {
        for (int j = offsets[i]; j < offsets[i + 1]; ++j)
        {
            int id_b = bond_ids[j];

            if (F[j * 4] != 0 || T[j * 4] != 0)
            {
                double val = cutoff[j * 2] * -attr_p[j * 2]
                             * (F[j * 4] + T[j * 4]);

                total_potential += val;

                potential[i] += val;
                potential[id_b] += val;

                bond_dir_force[j] += (F[j * 4] + T[j * 4]) * (cutoff[j * 2 + 1]
                                                              * -attr_p[j * 2] + cutoff[j * 2] * -attr_p[j * 2 + 1]);

                // dN_ti, dN_tj
                bond_dir_force[j] -= cutoff[j * 2] * -attr_p[j * 2]
                                     * (F[j * 4 + 1] + F[j * 4 + 2] + T[j * 4 + 1] + T[j * 4 + 2])
                                     * cutoff[j * 2 + 1];
                force_on_cutoff[i] += cutoff[j * 2] * -attr_p[j * 2]
                                      * (F[j * 4 + 1] + T[j * 4 + 1]);
                force_on_cutoff[id_b] += cutoff[j * 2] * -attr_p[j * 2]
                                         * (F[j * 4 + 2] + T[j * 4 + 2]);

                // dN_conj
                // dNca
                double mult = cutoff[j * 2] * -attr_p[j * 2]
                              * (F[j * 4 + 3] + T[j * 4 + 3]);

                for (int k = offsets[i]; k < offsets[i + 1]; ++k)
                {
                    int id_c = bond_ids[k];
                    if (cutoff[k * 2] == 0 ||
                        k == j ||
                        types[id_c] != atom_type::CARBON)
                        continue;

                    double Fc = F_conj(Nb_vars[k * 7 + 1]),
                        dFc = dF_conj(Nb_vars[k * 7 + 1]);
                    bond_dir_force[k] += mult * Nb_vars[j * 7 + 5]
                                         * cutoff[k * 2 + 1] * (Fc - cutoff[k * 2] * dFc);
                    force_on_cutoff[id_c] += mult * Nb_vars[j * 7 + 5]
                                             * cutoff[k * 2] * dFc;
                }

                // dNcb
                for (int k = offsets[id_b]; k < offsets[id_b + 1]; ++k)
                {
                    int id_c = bond_ids[k];
                    if (cutoff[k * 2] == 0 ||
                        k == sister_ids[j] ||
                        types[id_c] != atom_type::CARBON)
                        continue;

                    double Fc = F_conj(Nb_vars[k * 7 + 1]),
                        dFc = dF_conj(Nb_vars[k * 7 + 1]);
                    bond_dir_force[k] += mult * Nb_vars[j * 7 + 6]
                                         * cutoff[k * 2 + 1] * (Fc - cutoff[k * 2] * dFc);
                    force_on_cutoff[id_c] += mult * Nb_vars[j * 7 + 6]
                                             * cutoff[k * 2] * dFc;
                }
            }
        }
    }
}

void airebo::system_control_t::update_distribute()
{
    // references
    const vector<double> &lengths = bond_control.get_bond_lengths(),
        &deltas = bond_control.get_bond_deltas();
    const vector<int> &offsets = bond_control.get_bond_offsets(),
        &sister_ids = bond_control.get_sister_bonds();

    // cutoff dependent forces + bond direction forces
    for (int i = 0; i < na; ++i)
    {
        if (force_on_cutoff[i] != 0 && offsets[i + 1] > offsets[i])
        {
            for (int j = offsets[i]; j < offsets[i + 1]; ++j)
                bond_dir_force[j] += force_on_cutoff[i] * cutoff[2 * j + 1];
        }
        for (int j = offsets[i]; j < offsets[i + 1]; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                bond_force[j * 3 + k] += bond_dir_force[j]
                                         * deltas[j * 3 + k] / lengths[j];
                if (sister_ids[j] < j)
                {
                    bond_force[j * 3 + k] -= bond_force[sister_ids[j] * 3 + k];
                    bond_force[sister_ids[j] * 3 + k] = -bond_force[j * 3 + k];
                }
            }
        }
    }

    for (int i = 0; i < na; ++i)
        for (int j = offsets[i]; !static_atom[i] && j < offsets[i + 1]; ++j)
            for (int k = 0; k < 3; ++k)
                force[i * 3 + k] += bond_force[j * 3 + k];
}

diff_fn_t system_control_t::get_diff_fn()
{
    return [this](const auto &pos) {
        this->set_position(pos);
        this->update();

        return make_pair(this->get_value(), this->get_gradient());
    };
}

//void system_control_t::add_hydrogen()
//{
//    // distance to place hydrogen atoms off of carbon, in Angstroms
//    constexpr double ch_bond_len = 1.09;
//
//    // refresh bond components/lengths/etc
//    bond_control.update(position);
//
//    // get the bond graph, which represents atoms by vertices and bonds by edges
//    auto graph = bond_control.get_graph();
//    auto &deltas = bond_control.get_bond_deltas();
//
//    for (int atom_id = 0; atom_id < na; ++atom_id)
//    {
//        if (types[atom_id] == atom_type::HYDROGEN)
//            continue;
//
//        int n_bonds = graph.degree(atom_id);
//        switch (n_bonds)
//        {
//        case 0:
//            // no bonds, place hydrogen at (len, 0, 0) for reproducibility
//
//        case 1:
//
//        case 2:
//
//        default:
//            break;
//        }
//    }
//
//}
//

// TODO: make add_hydrogen potential agnostic
 void airebo::system_control_t::add_hydrogen()
 {
     int num_total = -1;
     int per_atom_max = -1;

     // update bond values
     update();

     // get references that will be needed later
     const vector<double> &lengths = bond_control.get_bond_lengths(),
         &deltas = bond_control.get_bond_deltas();
     const vector<int> &bond_ids = bond_control.get_bond_ids(),
         &offsets = bond_control.get_bond_offsets();

     // rank atoms by number of bonds
     vector<pair<double, int> > cvec;
     for (int i = 0; i < na; ++i)
         cvec.push_back(pair<double, int>(Na_vars[i * 4], i));
     sort(cvec.begin(), cvec.end());

     // add hydrogen
     int h_left = num_total;

     for (int i = 0; i < na; ++i)
     {
         int id = cvec[i].second;
         if (types[id] == atom_type::HYDROGEN)
             continue;

         if (Na_vars[id * 4] >= 3 && h_left <= 0)
             break;

         int c_hyd = 0;
         vector<pair<double, int> > bonds;
         for (int j = offsets[id]; j < offsets[id + 1]; ++j)
         {
             bonds.push_back(pair<double, int>(cutoff[j], j));
             if (types[bond_ids[j]] == atom_type::HYDROGEN && cutoff[j] >= 1)
                 ++c_hyd;
         }

         if (per_atom_max > 0 && c_hyd >= per_atom_max)
             continue;

         int nbonds = min((int)bonds.size(), 2);
         sort(bonds.begin(), bonds.end());

         vector<int> added_ids;
         for (int j = 0; j < nbonds; ++j)
             added_ids.push_back(bonds[j].second);

         // add hydrogen
         if (nbonds == 2)
         {
             int b_id_a = added_ids[0],
                 b_id_b = added_ids[1];

             double avg[3], len = 0;
             for (int j = 0; j < 3; ++j)
             {
                 avg[j] = deltas[b_id_a * 3 + j] / lengths[b_id_a]
                        + deltas[b_id_b * 3 + j] / lengths[b_id_b];
                 len += avg[j] * avg[j];
             }

             if (len == 0)
             {
                 auto avg_v = unit_normal_to(vec3_t(&deltas[b_id_a * 3]));
                 for (int k = 0; k < 3; ++k)
                     avg[k] = avg_v[k];

                 len = 1;
             }

             for (int j = 0; j < 3; ++j)
                 position.push_back(position[id * 3 + j] - 1.09 * avg[j]
                     / sqrt(len));
             types.push_back(atom_type::HYDROGEN);
             --h_left;
         }
         else if (nbonds == 1)
         {
             vec3_t bond_ref(&deltas[added_ids[0] * 3]);
             vec3_t axis = unit_normal_to(bond_ref);

             sp2::mat3x3_t rot;
             rot = util::gen_rotation(axis, 2 * M_PI / 3);
             vec3_t bond_a = (rot * bond_ref).unit_vector() * 1.09;

             rot = util::gen_rotation(axis, 4 * M_PI / 3);
             vec3_t bond_b = (rot * bond_ref).unit_vector() * 1.09;

             for (int j = 0; j < 3; ++j)
                 position.push_back(position[id * 3 + j] + 1.09 * bond_a[j]);
             types.push_back(atom_type::HYDROGEN);
             --h_left;

             if (per_atom_max > 1 && c_hyd == 0)
             {
                 for (int j = 0; j < 3; ++j)
                     position.push_back(position[id * 3 + j] + 1.09 * bond_b[j]);
                 types.push_back(atom_type::HYDROGEN);
                 --h_left;
             }
         }

         if (h_left == 0)
             break;
     }

     // set number of bonds and atoms
     na = types.size();
     ref_pos = position;

     // update
     update();
 }
