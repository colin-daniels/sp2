#include "common/neighbor/cell_array_t.hpp"

#include <iostream>
#include <cmath>

using namespace std;
using namespace sp2::fbc;

cell_array_t::cell_array_t() : n_cells(0), n_active(0),
    n_periodic(0), cell_lifetime(8),
    low_bound{}, high_bound{}, range{} {}

cell_array_t::~cell_array_t()
{
    for (auto ptr : cells_full)
        delete ptr;
}

void cell_array_t::init(int input_low[3], int input_high[3], int num_periodic)
{
    n_periodic = num_periodic;
    resize_array(input_low, input_high);
}

void cell_array_t::resize_array(int input_low[3], int input_high[3])
{
    int new_ncells = 1;
    bool need_resize = false;

    int input_range[3];
    for (int i = 0; i < 3; ++i)
    {
        new_ncells *= (input_high[i] - input_low[i]);
        input_range[i] = input_high[i] - input_low[i];

        if (i < n_periodic)
        {
            if (input_high[i] != high_bound[i] || input_low[i] != low_bound[i])
                need_resize = true;
        }
        else
        {
            // check if there were actual changes in the array dimensions, use a larger tolerance for shrinking the array
            if (input_high[i] - high_bound[i] < -2 || input_high[i] - high_bound[i] > 0 ||
                input_low[i] - low_bound[i] < 0 || input_low[i] - low_bound[i] > 2)
                need_resize = true;
        }
    }

    if (!need_resize)
        return;

    // make a new cell vector
    vector<periodic_cell_t*> new_cells_full(new_ncells, nullptr);

    // copy over existing cells
    for (int x = low_bound[0]; x < high_bound[0]; ++x)
    {
        for (int y = low_bound[1]; y < high_bound[1]; ++y)
        {
            for (int z = low_bound[2]; z < high_bound[2]; ++z)
            {
                bool in_new = true;
                if (z < input_low[2] || z >= input_high[2] ||
                    y < input_low[1] || y >= input_high[1] ||
                    x < input_low[0] || x >= input_high[0])
                    in_new = false;

                int id = index(x, y, z),
                    new_id = (x - input_low[0]) * input_range[1] * input_range[2]
                             + (y - input_low[1]) * input_range[2]
                             + (z - input_low[2]);

                // if the cell is in the new vector, copy it over
                if (in_new && cells_full[id] != nullptr)
                {
                    cells_full[id]->id = new_id;
                    for (unsigned int j = 0; j < cells_full[id]->atom_ids.size(); ++j)
                        last_indices[cells_full[id]->atom_ids[j]] = new_id;

                    new_cells_full[new_id] = cells_full[id];
                }
                else if (cells_full[id] != nullptr)
                {
                    // otherwise, set its atom indices to -1
                    for (unsigned int j = 0; j < cells_full[id]->atom_ids.size(); ++j)
                        last_indices[cells_full[id]->atom_ids[j]] = -1;

                    // and delete it
                    delete cells_full[id];
                    --n_active;
                }
            }
        }
    }

    // copy the new information
    n_cells = new_ncells;
    cells_full = new_cells_full;

    for (int i = 0; i < 3; ++i)
    {
        low_bound[i] = input_low[i];
        high_bound[i] = input_high[i];
        range[i] = input_range[i];
    }

    // relink all cells
    unlink_all();
    link_all();
}

void cell_array_t::unlink_all()
{
    for (int i = 0; i < n_cells; ++i)
        if (cells_full[i] != nullptr)
            cells_full[i]->unlink_all();
}

void cell_array_t::link_all()
{
    for (int x = low_bound[0]; x < high_bound[0]; ++x)
        for (int y = low_bound[1]; y < high_bound[1]; ++y)
            for (int z = low_bound[2]; z < high_bound[2]; ++z)
                link_cell(x, y, z);
}

void cell_array_t::link_cell(int id)
{
    int z = id % range[2],
        y = ((id - z) / range[2]) % range[1],
        x = ((id - z) / range[2] - y) / range[1];
    link_cell(x + low_bound[0], y + low_bound[1], z + low_bound[2]);
}

void cell_array_t::link_cell(int x, int y, int z)
{
    int cid = index(x, y, z);
    if (cells_full[cid] == nullptr)
        return;

    for (int dx = -1; dx < 2; ++dx)
    {
        for (int dy = -1; dy < 2; ++dy)
        {
            for (int dz = -1; dz < 2; ++dz)
            {
                int ix = ((x + dx) - low_bound[0] + range[0]) % range[0] + low_bound[0],
                    iy = ((y + dy) - low_bound[1] + range[1]) % range[1] + low_bound[1],
                    iz = ((z + dz) - low_bound[2] + range[2]) % range[2] + low_bound[2],
                    id = index(ix, iy, iz);

                int n_id = (dx + 1) * 9 + (dy + 1) * 3 + (dz + 1);

                if ((ix != x + dx && n_periodic < 1) ||
                    (iy != y + dy && n_periodic < 2) ||
                    (iz != z + dz && n_periodic < 3))
                    continue;

                // link
                if (n_id < 13 && cells_full[id] != nullptr)
                    cells_full[cid]->link(cells_full[id], n_id);
                else if (n_id > 13 && cells_full[id] != nullptr)
                    cells_full[id]->link(cells_full[cid], 26 - n_id);
            }
        }
    }
}

int cell_array_t::index(int x, int y, int z)
{
    return ((x - low_bound[0]) * range[1] + (y - low_bound[1])) * range[2] + (z - low_bound[2]);
}

void cell_array_t::update_cells(const vector<double> &lattice_positions, vector<double> &delta)
{
    // na is a max() to deal with adding/removing atoms from the system
    const int na = max<int>(lattice_positions.size() / 3, last_indices.size());
    vector<int> indices(na, -1);
    last_indices.resize(na, -1);
    last_icx.resize(na, -1);

    // get the actual indices
    get_indices(lattice_positions, indices);

    // clear cell check vectors and delete old cells
    for (int i = 0; i < n_cells; ++i)
    {
        if (cells_full[i] != nullptr)
        {
            cells_full[i]->update_clear();
            if (cells_full[i]->time_empty >= cell_lifetime)
            {
                delete cells_full[i];
                cells_full[i] = nullptr;
                --n_active;
            }
        }
    }

    // check for changes, add/remove atoms from cells
    for (int i = 0; i < na; ++i)
    {
        if (indices[i] != last_indices[i])
        {
            // remove the atom from its old cell, if it exists
            if (last_indices[i] < -1 || last_indices[i] >= n_cells || !std::isfinite(last_indices[i]))
                throw domain_error("atom index outside boundaries, last index: " + to_string(last_indices[i]));
            else if (last_indices[i] != -1 && cells_full[last_indices[i]] != nullptr)
                cells_full[last_indices[i]]->remove_atom(i, last_icx);

            // if an atom wasn't just removed from the system
            if (indices[i] < -1 || indices[i] >= static_cast<int>(cells_full.size()) || !std::isfinite(indices[i]))
                throw domain_error("atom index outside boundaries, index: " + to_string(indices[i]));
            else if (indices[i] != -1)
            {
                // if there is no cell to place the atom in, construct one
                if (cells_full[indices[i]] == nullptr)
                {
                    cells_full[indices[i]] = new periodic_cell_t(indices[i]);
                    ++n_active;

                    // link it
                    link_cell(indices[i]);
                }

                // add the atom to the cell
                cells_full[indices[i]]->add_atom(i, last_icx);
            }

            // update index
            last_indices[i] = indices[i];
        }

        if (i < static_cast<int>(delta.size()) && delta[i] == 1)
            cells_full[indices[i]]->add_check(i);
    }

    // get comparison vectors for all cells
    for (int i = 0; i < n_cells; ++i)
        if (cells_full[i] != nullptr)
            cells_full[i]->get_comp();

    // resize last indices down to the current number of atoms in the system
    last_indices.resize(lattice_positions.size() / 3);
    last_icx.resize(lattice_positions.size() / 3);
}


void cell_array_t::get_indices(const vector<double> &lattice_positions, vector<int> &indices)
{
    const int na = lattice_positions.size() / 3;

    // get indices
    int mult[3] = {range[2] * range[1], range[2], 1};
    for (int i = 0; i < na; ++i)
    {
        indices[i] = 0;
        for (int j = 0; j < 3; ++j)
            indices[i] += mult[j] * static_cast<int>(
                lattice_positions[i * 3 + j] - low_bound[j]);

        if (indices[i] < -1 || indices[i] > n_cells)
        {
            cout << i << ": out of bounds index: " << indices[i] << endl;
            for (int j = 0; j < 3; ++j)
            {
                cout << "(" << range[j] << ", " << low_bound[j] << ") ";
            }
            cout << endl;
        }
    }
}
