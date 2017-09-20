#include "common/structure_t.hpp"
#include "common/util/mpi.hpp"

#include "deps/json.hpp"
#include <common/math/rotations.hpp>
#include <common/math/vec3_util.hpp>

using namespace std;
using namespace sp2;
namespace mpi = boost::mpi;

structure_t::structure_t(const double lattice_in[3][3],
    const std::vector<atom_type> &types_in,
    const std::vector<double> &positions_in) :
        structure_t(lattice_in, types_in, sp2::dtov3(positions_in)) {}

structure_t::structure_t(const double lattice_in[3][3],
    const std::vector<atom_type> &types_in,
    const std::vector<vec3_t> &positions_in) :
        types(types_in),
        positions(positions_in)
{
    copy_n(lattice_in[0], 9, lattice[0]);
}

/// serialize object to Json::Value
bool structure_t::serialize(Json::Value &output) const
{
    if (!space_group.empty())
        output["space_group"] = space_group;

    output["lattice"] = serialize_lattice();

    Json::Value structure(Json::arrayValue);
    for (Json::Value::ArrayIndex i = 0; i < types.size(); ++i)
    {
        structure.append(Json::Value());
        structure[i].append(enum_to_str<atom_type>(types[i]));
        for (int j = 0; j < 3; ++j)
            structure[i].append(positions[i][j]);
    }

    output["positions"] = structure;
    return true;
}

/// deserialize object from Json::Value
bool structure_t::deserialize(const Json::Value &input)
{
    if (!input || !deserialize_lattice(input["lattice"]))
        return false;

    space_group = input.get("space_group", "").asString();

    if (!input["positions"] || !input["positions"].isArray())
        return false;

    types.clear();
    positions.clear();
    for (const Json::Value &atom_info : input["positions"])
    {
        if (!atom_info.isArray()     ||
            atom_info.size() != 4    ||
            !atom_info[0].isString())
            return false;

        string type = atom_info[0].asString();
        types.push_back(enum_from_str<atom_type>(type));

        sp2::vec3_t pos;
        for (int j = 1; j < 4; ++j)
        {
            if (!atom_info[j].isNumeric())
                return false;

            pos[j] = atom_info[j].asDouble();
        }

        positions.push_back(pos);
    }

    return true;
}

Json::Value structure_t::serialize_lattice() const
{
    Json::Value output;
    for (int i = 0; i < 3; ++i)
    {
        output.append(Json::Value());
        for (int j = 0; j < 3; ++j)
            output[i].append(lattice[i][j]);
    }
    return output;
}

bool structure_t::deserialize_lattice(const Json::Value &input)
{
    if (!input.isArray() || input.size() != 3)
        return false;

    for (int i = 0; i < 3; ++i)
    {
        if (!input[i].isArray() || input[i].size() != 3)
            return false;
        for (int j = 0; j < 3; ++j)
            lattice[i][j] = input[i].get(j, 0).asDouble();
    }
    return true;
}


void structure_t::bcast(const boost::mpi::communicator &comm, int root)
{
    mpi::broadcast(comm, lattice, root);
    mpi::broadcast(comm, positions, root);
    mpi::broadcast(comm, types, root);
}

void structure_t::rotate(vec3_t axis, double theta)
{
    transform(util::gen_rotation(axis, theta));
}

void structure_t::transform(const mat3x3_t &transformation)
{
    for (auto &p : positions)
        p = transformation * p;
}

// multiply many row vectors on the right by a matrix,
//  reusing allocated memory by 'out'
vector<vec3_t> right_multiply(const vector<vec3_t> &in,
    mat3x3_t mat, vector<vec3_t> out_buf)
{
    auto transpose = mat.transposed();

    auto out = move(out_buf);
    out.clear();
    out.reserve(in.size());
    for (auto vec : in)
        out.push_back(transpose * vec);
    return out;
}

// multiply many row vectors on the right by a matrix in-place
void right_multiply_inplace(vector<vec3_t> &vecs, mat3x3_t mat)
{
    auto transpose = mat.transposed();
    for (auto &vec : vecs)
        vec = transpose * vec;
}

std::vector<vec3_t> structure_t::fractional_positions() const
{
    auto inv = mat3x3_t(lattice).inverse();
    return right_multiply(positions, inv, {});
}

std::vector<vec3_t> structure_t::reduced_fractional_positions() const
{
    auto out = fractional_positions();
    for (auto &vec : out)
    {
        vec -= floor(vec); // into interval [0,1]; consider x = -1e-20
        vec -= floor(vec); // into interval [0,1)
    }
    return out;
}

void structure_t::set_fractional_positions(const vector<vec3_t> &fracs)
{
    positions = right_multiply(fracs, mat3x3_t{lattice}, move(positions));
}

void structure_t::reduce_positions()
{
    set_fractional_positions(reduced_fractional_positions());
}

vec3_t structure_t::get_lengths() const
{
    return {
        vec3_t(lattice[0]).mag(),
        vec3_t(lattice[1]).mag(),
        vec3_t(lattice[2]).mag()
    };
}

// This is what all of the (re)scale functions ultimately boil down to;
// modifying the lattice while preserving fractional positions.
void _set_lattice_constant_fracs(structure_t &structure, mat3x3_t new_lattice)
{
    // For constant frac:   x' . L'^{-1} == x . L^{-1}
    //                      x' == x . L^{-1} . L'
    auto latt_inv = mat3x3_t(structure.lattice).inverse();

    right_multiply_inplace(structure.positions, latt_inv * new_lattice);

    for (int i = 0; i < 3; i++)
        for (int k = 0; k < 3; k++)
            structure.lattice[i][k] = new_lattice[i][k];
}

void structure_t::scale(vec3_t factors)
{
    // OPTIMIZATION HAZARD:
    //
    // This does NOT boil down to multiplication of each column in positions
    //  by a constant, though such would appear to work for orthogonal cells.
    //
    // Cartesian positions in general transform as follows:
    // (x: position row vector, L: current lattice, s: scale factors)
    //
    //      x' == x . L^{-1} . diag(s) . L
    //
    scale({{factors[0], 0.0,        0.0},
           {0.0,        factors[1], 0.0},
           {0.0,        0.0,        factors[2]}});
}

void structure_t::scale(mat3x3_t factor)
{
    _set_lattice_constant_fracs(*this, factor * mat3x3_t(lattice));
}

void structure_t::rescale(vec3_t new_lengths)
{
    auto cur_lengths = get_lengths();
    return scale({
        new_lengths[0] / cur_lengths[0],
        new_lengths[1] / cur_lengths[1],
        new_lengths[2] / cur_lengths[2]
    });
}

void structure_t::rescale(mat3x3_t new_lattice)
{
    _set_lattice_constant_fracs(*this, new_lattice);
}

void sp2::sort_structure_types(structure_t &structure)
{
    vector<pair<atom_type, int>> atoms;
    for (std::size_t i = 0; i < structure.types.size(); ++i)
        atoms.emplace_back(structure.types[i], i);

    sort(atoms.begin(), atoms.end());
    sort(structure.types.begin(), structure.types.end());

    vector<vec3_t> sorted_pos;
    for (std::size_t i = 0; i < atoms.size(); ++i)
        sorted_pos.push_back(structure.positions[atoms[i].second]);

    structure.positions = sorted_pos;
}
