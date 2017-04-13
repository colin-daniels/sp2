#include "common/neighbor/bond_control_t.hpp"
#include "common/neighbor/utility_functions.hpp"
#include "common/math/blas.hpp"

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <common/math/vec3_t.hpp>

using namespace std;
using namespace sp2::fbc;

sp2::graph::ud_graph_t get_bond_graph(std::vector<double> positions,
    double lattice[3][3], double bond_max)
{
    bond_control_t bond_control;
    bond_control.init(lattice, bond_max, 0);

    bond_control.update(positions);
    return bond_control.get_graph();
}

void bond_control_t::init(const double input_lattice[3][3],
    double bond_max, double recalc_threshold)
{
    bool same = (cell_min     == bond_max + recalc_threshold) &&
                (cell_padding == recalc_threshold);

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            same = same && (input_lattice[i][j] == orig_lattice[i][j]);

    // dont do anything if nothing's changed
    if (same)
        return;

    cell_min = bond_max + recalc_threshold;
    cell_padding = recalc_threshold;

    set_lattice(input_lattice);
}

// set the lattice vectors for the system
void bond_control_t::set_lattice(const double input_lattice[3][3])
{
    // copy the input lattice
    copy_n(input_lattice[0], 9, orig_lattice[0]);

    // transpose to column major
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j)
            std::swap(orig_lattice[i][j], orig_lattice[j][i]);

    // process it and record the rotation applied to the lattice as well as the number of periodic directions
    n_periodic = process_lattice(orig_lattice, lattice, rotation);

    // scale the lattice vectors to cell widths
    int input_high[3], input_low[3];
    for (int i = 0; i < 3; ++i)
    {
        input_high[i] = cell_array.high_bound[i];
        input_low[i] = cell_array.low_bound[i];
    }

    for (int i = 0; i < 3; ++i)
    {
        // if its a periodic direction, ensure that the cells divide evenly
        double scaling = lattice[i][i] / cell_min;
        if (i < n_periodic)
        {
            input_low[i] = 0;
            input_high[i] = static_cast<int>(floor(scaling));
        }
        else
        {
            // scale for non-periodic directions
            for (int j = 0; j < 3; ++j)
                lattice[j][i] /= scaling;
        }

        if (input_low[i] == input_high[i])
            input_high[i] = input_low[i] + 1;
    }

    // get the inverse lattice matrix
    invert_3x3(lattice, inv_lattice);

    // and the transformation lattice along with its inverse
    //     transformation = rot * inv_lattice
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            transformation[i][j] = 0;
            for (int k = 0; k < 3; ++k)
                transformation[i][j] += inv_lattice[i][k] * rotation[k][j];
        }
    }
    invert_3x3(transformation, inv_transformation);

    // initialize the cell array
    cell_array.init(input_low, input_high, n_periodic);
}

void bond_control_t::set_bonds(const vector<int> &input_offsets, const vector<int> &input_bond_ids)
{
    for (unsigned int i = 0; i < full_bond_list.size(); ++i)
        for (unsigned int j = 0; j < full_bond_list[i].size(); ++j)
            full_bond_list[i][j] = false;

    offsets = input_offsets;
    bond_ids = input_bond_ids;

    for (int i = 0; i < static_cast<int>(offsets.size()) - 1; ++i)
        for (int j = offsets[i]; j < offsets[i + 1]; ++j)
            full_bond_list[i][bond_ids[j]] = true;

    nb = bond_ids.size();
}

void bond_control_t::resize_members(int num_atoms)
{
    if (num_atoms != na)
    {
        na = num_atoms;
        last_positions.resize(na * 3, 1e99);
        for (int i = 0; i < min(na, static_cast<int>(full_bond_list.size())); ++i)
            full_bond_list[i].resize(na, false);
        full_bond_list.resize(na, vector<bool>(na, false));

        offsets.resize(na + 1, 0);
    }
}

void bond_control_t::mark_movement(const vector<double> &positions, vector<double> &delta)
{
    delta = positions;
    for (auto i = 0u; i < positions.size(); ++i)
        delta[i] -= last_positions[i];

    // check position changes
    for (int i = 0; i < na; ++i)
    {
        double dist = delta[i*3] * delta[i*3]
                      + delta[i*3 + 1] * delta[i*3 + 1]
                      + delta[i*3 + 2] * delta[i*3 + 2];

        // mark them
        if (dist > (cell_padding * cell_padding / 4.0))
            delta[i] = 1;
        else if (!std::isfinite(delta[i]))
            throw domain_error("Error, non-finite atom position.");
        else
            delta[i] = 0;
    }
    delta.resize(na);
}

void bond_control_t::wrap_positions(vector<double> &positions)
{
    // transform into lattice coordinates
    vector<double> lattice_positions(na * 3, 0);
    loop_dgemm<3, 3>(positions.data(), na, transformation[0], lattice_positions.data());

    // wrap in the periodic directions
    for (int i = 0; i < n_periodic; ++i)
        for (int j = 0; j < na; ++j)
            lattice_positions[j*3 + i] -= floor(lattice_positions[j*3 + i]);

    // transform back to normal coordinates
    loop_dgemm<3, 3>(lattice_positions.data(), na, inv_transformation[0], positions.data());
}

void bond_control_t::calc_delta(const double *pos_a, const double *pos_b, double *output) const
{
    double delta[3] = {};
    for (int i = 0; i < 3; ++i)
        delta[i] = pos_b[i] - pos_a[i];

    // transform positions into lattice coordinates
    double lattice_delta[3] = {};
    loop_dgemv<3, 3>(transformation[0], delta, lattice_delta);

    // wrap
    for (int i = 0; i < n_periodic; ++i)
        lattice_delta[i] -= round(lattice_delta[i]);

    // transform delta back into normal coordinates
    loop_dgemv<3, 3>(inv_transformation[0], lattice_delta, output);
}

void bond_control_t::calc_delta_all(std::vector<double> &deltas, std::vector<double> &working)
{
    // get lattice delta values
    deltas.resize(3 * ((na * na) / 2));

    for (size_t i = 0, tpos = 0; i < last_lattice_positions.size(); i += 3)
    {
        for (size_t j = i + 3; j < last_lattice_positions.size(); j += 3, tpos += 3)
        {
            for (int k = 0; k < 3; ++k)
            {
                deltas[tpos + k] = last_lattice_positions[j + k] - last_lattice_positions[i + k];
                if (k < n_periodic)
                    deltas[tpos + k] -= round(deltas[tpos + k]);
            }
        }
    }

    // transform deltas to normal coordinates
    transform_from_lattice(deltas, working);
    swap(deltas, working);
}

std::vector<double> bond_control_t::calc_delta_graph(const graph::ud_graph_t &input_graph,
    const std::vector<double> &input_positions)
{
    // transform atom positions to lattice positions
    vector<double> lattice_positions;
    transform_to_lattice(input_positions, lattice_positions);

    // get deltas
    vector<double> lattice_deltas(input_graph.n_edges() * 3);
    for (auto edge : input_graph.edges())
    {
        double *delta = &lattice_deltas[edge.id * 3];
        for (int k = 0; k < 3; ++k)
            delta[k] = lattice_positions[edge.b * 3 + k]
                       - lattice_positions[edge.a * 3 + k];

        // wrap periodic coordinates
        for (int k = 0; k < n_periodic; ++k)
            delta[k] -= round(delta[k]);
    }

    // transform deltas back into normal coordinates
    vector<double> output_deltas;
    transform_from_lattice(lattice_deltas, output_deltas);

    return output_deltas;
}

double bond_control_t::calc_distance(const double *pos_a, const double *pos_b) const
{
    // get the delta values
    double delta[3];
    calc_delta(pos_a, pos_b, delta);

    // return length
    return sqrt(delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]);
}

void bond_control_t::transform_to_lattice(const vector<double> &positions, vector<double> &lattice_positions) const
{
    lattice_positions.resize(positions.size());
    loop_dgemm<3, 3>(positions.data(), positions.size() / 3, transformation[0], lattice_positions.data());
}

void bond_control_t::transform_from_lattice(const vector<double> &lattice_positions, vector<double> &positions) const
{
    positions.resize(lattice_positions.size());
    loop_dgemm<3, 3>(lattice_positions.data(), lattice_positions.size() / 3,
        inv_transformation[0], positions.data());
}

void bond_control_t::update(const vector<double> &positions)
{
    // resize members if needed
    resize_members(positions.size() / 3);
    if (na <= 1)
        return;

    // get the indices of atoms that moved a lot
    vector<double> delta;
    mark_movement(positions, delta);

    vector<double> lattice_positions;
    transform_to_lattice(positions, lattice_positions);

    int range[3];
    cell_array.get_range(range);

    // wrap periodic coordinates
    for (int i = 0; i < n_periodic; ++i)
    {
        for (int j = 0; j < na; ++j)
        {
            auto &pos = lattice_positions[j * 3 + i];
            pos = min(pos - floor(pos), nexttoward(1.0, 0));
        }
    }

    // record positions
    last_lattice_positions = lattice_positions;

    // special update for small systems
    if (na < 50)
    {
        update_small();
        return;
    }

    // scale periodic coordinates
    double drange[3] = {1.0 * range[0], 1.0 * range[1], 1.0 * range[2]};
    for (auto i = 0u; i < lattice_positions.size(); i += 3)
        for (auto j = 0; j < n_periodic; ++j)
            lattice_positions[i + j] *= drange[j];

    // get non-periodic min/max values for cell array resizing
    int min_index[3] = {},
        max_index[3] = {};
    for (int i = 0; i < n_periodic; ++i)
        max_index[i] = range[i];
    for (int i = n_periodic; i < 3; ++i)
    {
        for (int j = 0; j < na; ++j)
        {
            if (j == 0)
            {
                min_index[i] = static_cast<int>(floor(lattice_positions[j*3 + i]));
                max_index[i] = static_cast<int>(ceil(lattice_positions[j*3 + i]));
            }
            else
            {
                min_index[i] = min(min_index[i], static_cast<int>(floor(lattice_positions[j*3 + i])));
                max_index[i] = max(max_index[i], static_cast<int>(ceil(lattice_positions[j*3 + i])));
            }
        }
        if (min_index[i] == max_index[i])
            ++max_index[i];
    }

    // resize if needed
    cell_array.resize_array(min_index, max_index);

    // update cells and get indices
    cell_array.update_cells(lattice_positions, delta);

    // get bonds
    vector<double> temp_deltas;
    vector<pair<int, int> > temp_bonds;
    get_bonds(lattice_positions, temp_deltas, temp_bonds);

    // process bonds
    process_bonds(temp_deltas, temp_bonds);

    for (int i = 0; i < na; ++i)
        if (delta[i] == 1)
            for (int j = 0; j < 3; ++j)
                last_positions[i * 3 + j] = positions[i * 3 + j];
}

void bond_control_t::process_bonds(vector<double> &temp_deltas, vector<pair<int, int> > &temp_bonds)
{
    // transform to regular coordinates
    vector<double> temp = temp_deltas;
    loop_dgemm<3, 3>(temp.data(), temp.size() / 3, inv_transformation[0], temp_deltas.data());

    vector<int> prc_offsets(na + 1, 0);
    vector<double> temp_lengths;
    for (unsigned int i = 0; i < temp_deltas.size(); i += 3)
    {
        temp_lengths.push_back(temp_deltas[i] * temp_deltas[i] + temp_deltas[i+1] * temp_deltas[i+1] + temp_deltas[i+2] * temp_deltas[i+2]);
        int id_a = temp_bonds[i/3].first,
            id_b = temp_bonds[i/3].second;
        if (id_a == id_b)
            continue;

        if (bond_lock
            || (temp_lengths.back() < cell_min * cell_min
                || (static_cast<int>(i) < (3 * nb / 2)
                    && temp_lengths.back() < (cell_min + cell_padding) * (cell_min + cell_padding))))
        {
            ++prc_offsets[id_a + 1];
            ++prc_offsets[id_b + 1];
        }
        else
        {
            full_bond_list[id_a][id_b] = false;
            full_bond_list[id_b][id_a] = false;
        }
    }

    for (int i = 0; i < na; ++i)
        prc_offsets[i + 1] += prc_offsets[i];

    int old_nb = nb;
    nb = prc_offsets.back();
    vector<int> prc_bonds(nb, -1),
        prc_sister_ids(nb);
    vector<double> prc_deltas(nb* 3),
        prc_lengths(nb);

    vector<int> last_input = prc_offsets;
    for (unsigned int i = 0; i < temp_bonds.size(); ++i)
    {
        if (temp_bonds[i].first == temp_bonds[i].second)
            continue;

        if (bond_lock || (temp_lengths[i] < cell_min * cell_min || (static_cast<int>(i) < old_nb / 2 && temp_lengths[i] < (cell_min + cell_padding) * (cell_min + cell_padding))))
        {
            int id_a = temp_bonds[i].first,
                id_b = temp_bonds[i].second;

            prc_bonds[last_input[id_a]] = id_b;
            prc_bonds[last_input[id_b]] = id_a;

            prc_sister_ids[last_input[id_a]] = last_input[id_b];
            prc_sister_ids[last_input[id_b]] = last_input[id_a];

            prc_lengths[last_input[id_a]] = sqrt(temp_lengths[i]);
            prc_lengths[last_input[id_b]] = sqrt(temp_lengths[i]);

            for (int k = 0; k < 3; ++k)
            {
                prc_deltas[last_input[id_a] * 3 + k] =  temp_deltas[i * 3 + k];
                prc_deltas[last_input[id_b] * 3 + k] = -temp_deltas[i * 3 + k];
            }

            ++last_input[id_a];
            ++last_input[id_b];
        }
    }

    swap(offsets, prc_offsets);
    swap(bond_ids, prc_bonds);
    swap(delta_list, prc_deltas);
    swap(length_list, prc_lengths);
    swap(sister_ids, prc_sister_ids);
}

void bond_control_t::get_bonds(vector<double> &lattice_positions, vector<double> &temp_deltas, vector<pair<int, int> > &temp_bonds)
{
    // scale back periodic coordinates
    int range[3] = {};
    cell_array.get_range(range);
    double drange[3] = {1.0 / range[0], 1.0 / range[1], 1.0 / range[2]};
    for (auto i = 0u; i < lattice_positions.size(); i += 3)
        for (auto j = 0; j < n_periodic; ++j)
            lattice_positions[i + j] *= drange[j];

    // exclude bonds via the full_bond_list array, first set to false
    for (unsigned int i = 0; i < excluded_bonds.size(); ++i)
    {
        if (excluded_bonds[i].first >= na || excluded_bonds[i].second >= na)
            continue;
        full_bond_list[excluded_bonds[i].first][excluded_bonds[i].second] = false;
        full_bond_list[excluded_bonds[i].second][excluded_bonds[i].first] = false;
    }

    // get existing bonds
    for (int id_a = 0; id_a < na; ++id_a)
    {
        for (int j = offsets[id_a]; j < offsets[id_a + 1]; ++j)
        {
            if (id_a > bond_ids[j] && full_bond_list[id_a][bond_ids[j]])
            {
                temp_bonds.push_back(pair<int, int>(id_a, bond_ids[j]));
                for (int k = 0; k < 3; ++k)
                {
                    temp_deltas.push_back(lattice_positions[bond_ids[j] * 3 + k] - lattice_positions[3 * id_a + k]);
                    if (k < n_periodic)
                        temp_deltas.back() -= round(temp_deltas.back());
                }
            }
        }
    }

    // finish excluding
    for (unsigned int i = 0; i < excluded_bonds.size(); ++i)
    {
        if (excluded_bonds[i].first >= na || excluded_bonds[i].second >= na)
            continue;
        full_bond_list[excluded_bonds[i].first][excluded_bonds[i].second] = true;
        full_bond_list[excluded_bonds[i].second][excluded_bonds[i].first] = true;
    }

    if (bond_lock)
        return;

    // new bonds
    for (int c = 0; c < cell_array.size(); ++c)
    {
        if (cell_array.cells_full[c] == nullptr)
            continue;
        periodic_cell_t &cell = *cell_array.cells_full[c];

        for (unsigned int i = 0; i < cell.atom_ids.size(); ++i)
        {
            int id_a = cell.atom_ids[i];

            // new bonds #1
            for (unsigned int j = 0; j < cell.to_check.size(); ++j)
            {
                if (cell.to_check[j] != id_a && !full_bond_list[id_a][cell.to_check[j]])
                {
                    temp_bonds.push_back(pair<int, int>(id_a, cell.to_check[j]));
                    for (int k = 0; k < 3; ++k)
                    {
                        temp_deltas.push_back(lattice_positions[cell.to_check[j] * 3 + k] - lattice_positions[3 * id_a + k]);
                        if (k < n_periodic)
                            temp_deltas.back() -= round(temp_deltas.back());
                    }
                    full_bond_list[id_a][cell.to_check[j]] = true;
                    full_bond_list[cell.to_check[j]][id_a] = true;
                }
            }

            // new bonds #2
            for (unsigned int j = 0; j < cell.comp.size(); ++j)
            {
                if (!full_bond_list[id_a][cell.comp[j]])
                {
                    temp_bonds.push_back(pair<int, int>(id_a, cell.comp[j]));
                    for (int k = 0; k < 3; ++k)
                    {
                        temp_deltas.push_back(lattice_positions[cell.comp[j] * 3 + k] - lattice_positions[3 * id_a + k]);
                        if (k < n_periodic)
                            temp_deltas.back() -= round(temp_deltas.back());
                    }
                    full_bond_list[id_a][cell.comp[j]] = true;
                    full_bond_list[cell.comp[j]][id_a] = true;
                }
            }
        }

        for (unsigned int i = 0; i < cell.to_check.size(); ++i)
        {
            int id_a = cell.to_check[i];

            // new bonds #3
            for (int j = 0; j < 13; ++j)
            {
                if (cell.neighbors[j] != nullptr)
                {
                    for (unsigned int l = 0; l < cell.neighbors[j]->atom_ids.size(); ++l)
                    {
                        int id_b = cell.neighbors[j]->atom_ids[l];
                        if (!full_bond_list[id_a][id_b])
                        {
                            temp_bonds.push_back(pair<int, int>(id_a, id_b));
                            for (int k = 0; k < 3; ++k)
                            {
                                temp_deltas.push_back(lattice_positions[id_b * 3 + k] - lattice_positions[3 * id_a + k]);
                                if (k < n_periodic)
                                    temp_deltas.back() -= round(temp_deltas.back());
                            }
                            full_bond_list[id_a][id_b] = true;
                            full_bond_list[id_b][id_a] = true;
                        }
                    }
                }
            }
        }
    }
}


void bond_control_t::set_excluded_bond(int id_a, int id_b)
{
    pair<int, int> bond(id_a, id_b);
    if (id_a < id_b)
        swap(bond.first, bond.second);
    excluded_bonds.push_back(bond);
}


void bond_control_t::remove_excluded_bond(int id_a, int id_b)
{
    pair<int, int> bond(id_a, id_b);
    if (id_a < id_b)
        swap(bond.first, bond.second);
    vector<pair<int, int> >::iterator it = find(excluded_bonds.begin(), excluded_bonds.end(), bond);
    if (it != excluded_bonds.end())
        excluded_bonds.erase(it);
}

void bond_control_t::clear_excluded_bonds()
{
    excluded_bonds.clear();
}

void bond_control_t::update_small()
{
    std::vector<double> prc_deltas,
        prc_lengths;

    std::vector<int> prc_bonds,
        prc_sister_ids,
        prc_offsets = {0};

    const double max_dist = (cell_min - cell_padding)
                            * (cell_min - cell_padding);

    const auto& llp = last_lattice_positions;
    auto get_delta = [&](int i, int j) {
        vec3_t delta = {
            llp[j * 3 + 0] - llp[i * 3 + 0],
            llp[j * 3 + 1] - llp[i * 3 + 1],
            llp[j * 3 + 2] - llp[i * 3 + 2]
        };

        // wrap for periodic directions
        for (int k = 0; k < n_periodic; ++k)
            delta[k] -= round(delta[k]);

        return delta;
    };

    auto update_prc = [&](const auto &lattice_delta, int i, int j, int &m) {
        auto delta = lattice_delta.mul_3x3(inv_transformation);
        auto len_sq = delta.mag_sq();
        if (len_sq > max_dist)
            return false;

        prc_bonds.push_back(j);
        prc_sister_ids.push_back(-1);

        for (auto v : delta)
            prc_deltas.push_back(v);

        prc_lengths.push_back(sqrt(len_sq));

        m += 1;
        if (j > i)
            return true;

        int nm = 0;
        for (int k = prc_offsets[j]; k < prc_offsets[j + 1]; ++k)
        {
            if (prc_bonds[k] == i && nm++ == (m - 1))
            {
                prc_sister_ids.back() = k;
                prc_sister_ids[k] = prc_bonds.size() - 1;
                break;
            }
        }

        return true;
    };

    for (size_t i = 0; i < llp.size() / 3; ++i)
    {
        for (size_t j = 0; j < llp.size() / 3; j++)
        {
            if (i == j)
                continue;

            int m = 0;
            auto lattice_delta = get_delta(i, j);
            if (!update_prc(lattice_delta, i, j, m))
                continue;

            auto modify = [&](vector<int> dim) {
                auto ldelta_mod = lattice_delta;
                for (int k = 0; k < 3; ++k)
                    if (dim[k])
                        ldelta_mod[k] += ldelta_mod[k] < 0 ? 1 : -1;

                return ldelta_mod;
            };

            if (n_periodic == 0)
                continue;

            update_prc(modify({1, 0, 0}), i, j, m);
            if (n_periodic == 1)
                continue;

            update_prc(modify({0, 1, 0}), i, j, m);
            update_prc(modify({1, 1, 0}), i, j, m);
            if (n_periodic == 2)
                continue;

            update_prc(modify({1, 1, 1}), i, j, m);
            update_prc(modify({1, 0, 1}), i, j, m);
            update_prc(modify({0, 1, 1}), i, j, m);
            update_prc(modify({0, 0, 1}), i, j, m);
        }
        prc_offsets.push_back(prc_bonds.size());
    }

    offsets = move(prc_offsets);
    bond_ids = move(prc_bonds);
    delta_list = move(prc_deltas);
    length_list = move(prc_lengths);
    sister_ids = move(prc_sister_ids);

    nb = bond_ids.size();
    for (std::size_t i = 0; i < sister_ids.size(); ++i)
        assert (sister_ids[sister_ids[i]] == static_cast<int>(i));
}
