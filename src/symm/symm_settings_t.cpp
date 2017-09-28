#include "symm/symm_settings_t.hpp"
#include "common/json/json.hpp"

using namespace std;
using namespace sp2;

void symm::symm_settings_t::serialize(Json::Value &output) const
{
    io::serialize_basic(output,
        "n_atoms", n_atoms,
        "space_group", space_group,
        "space_group_file", space_group_file,
        "n_connected", n_connected,
        "output_filename", output_filename,
        "use_cg", use_cg,
        "unit_cell_range", unit_cell_range,
        "bond_range", bond_range,
        "pso", pso_set,
        "acgsd", acgsd_set
    );
}

bool symm::symm_settings_t::deserialize(const Json::Value &input)
{
    return io::deserialize_basic(input,
        "n_atoms", n_atoms,
        "space_group", space_group,
        "space_group_file", space_group_file,
        "n_connected", n_connected,
        "output_filename", output_filename,
        "use_cg", use_cg,
        "unit_cell_range", unit_cell_range,
        "bond_range", bond_range,
        "pso", pso_set,
        "acgsd", acgsd_set
    );
}
