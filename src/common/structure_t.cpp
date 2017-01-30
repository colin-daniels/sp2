#include "common/structure_t.hpp"
#include "common/util/mpi.hpp"

#include <json/json.h>

using namespace std;
using namespace sp2;
namespace mpi = boost::mpi;

structure_t::structure_t(const double lattice_in[3][3],
    const std::vector<atom_type> &types_in,
    const std::vector<double> &positions_in) :
    types(types_in), positions(positions_in)
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
            structure[i].append(positions[i * 3 + j]);
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

        for (int j = 1; j < 4; ++j)
        {
            if (!atom_info[j].isNumeric())
                return false;
            positions.push_back(atom_info[j].asDouble());
        }
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

void sp2::sort_structure_types(structure_t &structure)
{
    vector<pair<atom_type, int>> atoms;
    for (int i = 0; i < structure.types.size(); ++i)
        atoms.emplace_back(structure.types[i], i);

    sort(atoms.begin(), atoms.end());
    sort(structure.types.begin(), structure.types.end());

    vector<double> sorted_pos;
    for (int i = 0; i < atoms.size(); ++i)
        for (int j = 0; j < 3; ++j)
            sorted_pos.push_back(structure.positions[atoms[i].second * 3 + j]);

    structure.positions = sorted_pos;
}
