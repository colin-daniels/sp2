#ifndef SP2_NEIGHBOR_CONTROL_T_HPP
#define SP2_NEIGHBOR_CONTROL_T_HPP

#include <vector>
#include <stdexcept>
#include <algorithm>

#include "common/vec3_t.hpp"
#include "common/graph/ud_graph_t.hpp"

namespace sp2 {
namespace nc {
/// single periodic cell class
struct periodic_cell_t
{


};

class cell_array_t
{
private:
    int low_bound[3],
        high_bound[3],
        range[3];

    int n_periodic;
    std::vector<int> indices;

    std::vector<periodic_cell_t*> cells;
public:

    void resize(const int input_low[3], const int input_high[3],
        int n_periodic_in)
    {
        const int input_range[3] = {
            input_high[0] - input_low[0],
            input_high[1] - input_low[1],
            input_high[2] - input_low[2]
        };

        const int input_size = input_range[0] * input_range[1] * input_range[2];

        std::vector<periodic_cell_t*> new_cells(input_size, nullptr);

    }

    // if a particular dimension is considered periodic
    bool dim_periodic(int dim)
    {
        return dim < n_periodic;
    }

    void update(const std::vector<vec3_t> &positions)
    {
        if (positions.empty())
            return;

        // check if a resize is needed, and if so, do it
        check_and_resize(positions);

        // get indices
        int mult[3] = {range[2] * range[1], range[1], 1};
        for (int i = 0; i < n_periodic; ++i)
            mult[i] *= range[i];

        std::vector<int> indices;
        indices.reserve(positions.size());

        for (const auto &v : positions)
        {
            int index = 0;
            for (int i = 0; i < n_periodic; ++i)
                index += mult[i] * static_cast<int>(v[i]);
            for (int i = n_periodic; i < 3; ++i)
                index += mult[i] * static_cast<int>(v[i] - low_bound[i]);

            indices.push_back(index);
        }
    }


private:
    void check_and_resize(const std::vector<vec3_t> &positions)
    {
        // assumes positions is not empty
        auto min_v = positions.front(),
            max_v = min_v;

        for (auto &v : positions)
        {
            min_v = min_elem(min_v, v);
            max_v = min_elem(max_v, v);
        }

        int input_low[3], input_high[3];
        std::copy_n(low_bound,  3, input_low);
        std::copy_n(high_bound, 3, input_high);

        // check for out of bounds positions in periodic directions
        for (int i = 0; i < n_periodic; ++i)
            if (min_v[i] < low_bound[i] ||
                max_v[i] >= high_bound[i])
                throw std::domain_error("vertex out of cell array bounds");

        // determine if we need to resize the cell array in
        // non-periodic directions
        bool need_resize = false;
        for (int i = n_periodic; i < 3; ++i)
        {
            if (min_v[i] < low_bound[i]     ||
                max_v[i] >= high_bound[i]   ||
                // for shrinking the cell array we use a larger tolerance
                // to avoid oscillations if one vertex goes back and forth
                // across a cell boundary
                min_v[i] > low_bound[i] + 1 ||
                max_v[i] < high_bound[i] - 2)
            {
                need_resize = true;
            }

            input_low[i]  = static_cast<int>(std::floor(min_v[i]));
            input_high[i] =  static_cast<int>(std::ceil(max_v[i]));
        }

        if (need_resize)
            resize(input_low, input_high, n_periodic);
    }
};

class neighbor_control_t
{
private:

    /// previous positions passed to update(), used to only
    /// update vertices that move a significant amount
    std::vector<vec3_t> last_positions;

    /// number of periodic directions
    int n_periodic;
    /// normal coords -> lattice coords
    double lat_transform[3][3],
    /// lattice coords -> normal coords
        inv_transform[3][3];

    /// graph describing the adjacency of the vertex positions
    /// last passed to the update() function
    graph::ud_graph_t graph;

    cell_array_t cell_array;

    double recalc_threshold;

public:
    neighbor_control_t() = default;

    // non-copyable
    neighbor_control_t(const neighbor_control_t&) = delete;
    neighbor_control_t& operator=(const neighbor_control_t&) = delete;

    void init(const double lattice[3][3], double max_dist,
        double recalc_threshold)
    {

    }

    void update(const std::vector<vec3_t> &positions)
    {
        const auto nv = positions.size();

        // check inputs and mark vertices that need recalculation
        std::vector<bool> to_update(nv, true);
        for (auto i = 0u; i < nv; ++i)
        {
            for (auto d : positions[i])
                if (!std::isfinite(d))
                    throw std::domain_error("Non-finite input position.");

            if (i > last_positions.size())
                continue;

            double dist = (positions[i] - last_positions[i]).mag_sq();
            if (dist < recalc_threshold * recalc_threshold)
                to_update[i] = false;
        }

        // convert positions to lattice coordinates
        auto lattice_positions = positions;
        to_lattice(lattice_positions);

        // wrap in periodic directions, use std::min() to
        // avoid (v[i] - floor(v[i])) == 1 for small negative v[i]
        for (auto &v : lattice_positions)
            for (auto i = 0; i < n_periodic; ++i)
                v[i] = std::min<double>(v[i] - std::floor(v[i]),
                    std::nexttoward(1.0, 0));

        // update cell array



        // record new positions
        last_positions.resize(nv);
        for (auto i = 0u; i < nv; ++i)
            if (to_update[i])
                last_positions[i] = positions[i];
    }

    /// normal coords -> lattice coords
    void to_lattice(std::vector<vec3_t> &positions) const
    {
        for (auto &v : positions)
            v = v.mul_3x3(lat_transform);
    }

    /// lattice coords -> normal coords
    void from_lattice(std::vector<vec3_t> &positions) const
    {
        for (auto &v : positions)
            v = v.mul_3x3(inv_transform);
    }

    graph::ud_graph_t get_graph() const {
        return graph;}

private:
    void process_lattice(const double lattice[3][3], double cutoff)
    {
        // copy the input lattice for manipulation
        std::array<vec3_t, 3> lat_vec;
        for (int i = 0; i < 3; ++i)
            lat_vec[i] = vec3_t(lattice[i]);

        // non-zero lattice vectors are periodic directions, put them
        // at the "top" of the array and record how many there are
        int n_non_zero = 0;
        for (int i = 0; i < 3; ++i)
            if (lat_vec[i].mag() > 0)
                std::swap(lat_vec[i], lat_vec[n_non_zero++]);

        // set any zero vectors to be orthogonal to the non-zero ones
        // note: switch fall-through is intended
        switch (n_non_zero)
        {
        case 0: lat_vec[0] = vec3_t(1, 0, 0);
        case 1: lat_vec[1] = unit_normal_to(lat_vec[0]);
        case 2: lat_vec[2] = cross(lat_vec[0], lat_vec[1]);
        case 3:
        default:
            break;
        }

        // scale non-periodic directions
        for (int i = n_non_zero; i < 3; ++i)
            lat_vec[i] = lat_vec[i].unit_vector() * cutoff;

        // get 3x3 inverse via cross products
        std::array<vec3_t, 3> inv_vec;
        for (int i = 0; i < 3; ++i)
            inv_vec[i] = cross(lat_vec[(i + 1) % 3], lat_vec[(i + 2) % 3]);

        // check lattice vectors for linear dependence (determinant)
        double det = dot(lat_vec[0], inv_vec[0]);
        if (det == 0)
            throw std::domain_error(
                "Error, linearly dependent lattice vectors.");
        else
        {
            // finish 3x3 inverse calculation
            for (auto &v : inv_vec)
                v /= det;
        }

        // transpose and record the transformation matrices
        n_periodic = n_non_zero;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                // lattice coords -> normal coords
                inv_transform[j][i] = lat_vec[i][j];

                // normal coords -> lattice coords
                lat_transform[j][i] = inv_vec[i][j];
            }
        }

        // set up cell array boundaries
        int low_bound[3]  = {0, 0, 0},
            high_bound[3] = {1, 1, 1};

        // get scaling for periodic directions
        // a = cutoff / (l1 * (l2 x l3) / || l2 x l3 ||)
        // where a * l1 is the min len of l1 so that a sphere of cutoff
        // radius can fit inside the unit cell in the l1 direction
        for (int i = 0; i < n_periodic; ++i)
            high_bound[i] = std::floor(cutoff / inv_vec[i].mag());

        cell_array.resize(low_bound, high_bound, n_periodic);
    }
};

} // namespace nc
} // namespace sp2

#endif // SP2_NEIGHBOR_CONTROL_T_HPP
