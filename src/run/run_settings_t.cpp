#include "common/json/json.hpp"
#include "run/run_settings_t.hpp"

using namespace std;
using namespace sp2;

bool load_structure(const Json::Value &config, structure_t &structure)
{
    // structural information can either be in the same file as the rest
    // of the configuration data (aka config.json) or can be specified as
    // being in a separate file using a field in config.json containing the
    // file name.

    // first check if a separate filename is specified
    string structure_file = config.get("structure_file", "").asString();
    if (!structure_file.empty())
    {
        if (!io::read_structure(structure_file, structure))
        {
            cerr << "Error, failed to read structure file." << endl;
            return false;
        }
    }
    else
    {
        // if no separate file was specified, check for data within the
        // configuration file itself
        if (!io::deserialize_field(config, structure, "structure"))
            return false;
    }
    return true;
}

bool run_settings_t::serialize(Json::Value &output) const
{
    output["run_type"] = "";
    output["structure_file"] = "";
    output["potential_type"] = "";

    io::serialize_basic(output,
        "add_hydrogen", add_hydrogen,
        "log_filename", log_filename,
        "n_threads", n_threads,
        "symm", symm_settings,
        "atac", atac_settings,
#ifdef SP2_ENABLE_LAMMPS
        "lammps", lammps_settings,
#endif // SP2_ENABLE_LAMMPS
#ifdef SP2_ENABLE_PHONOPY
        "phonopy", phonopy_settings,
#endif // SP2_ENABLE_PHONOPY
        "relax", relaxation_settings
    );

    return true;
}

bool run_settings_t::deserialize(const Json::Value &input)
{
    mode = io::deserialize_enum<run_type>(input, "run_type", {
        {"atac",    run_type::ATAC},  // acceptable inputs
        {"relax",   run_type::RELAX}, //
        {"symm",    run_type::SYMM},  //
#ifdef SP2_ENABLE_PHONOPY
        {"phonopy", run_type::PHONOPY}
#endif // SP2_ENABLE_PHONOPY
    }, run_type::NONE, true);

    potential = io::deserialize_enum<potential_type>(input, "potential_type", {
        {"rebo",   potential_type::REBO},  // acceptable inputs
        {"lammps", potential_type::LAMMPS} //
    }, potential_type::NONE, true);

    if (mode == run_type::NONE)
        return false;

    if (potential == potential_type::NONE)
        return false;

    if (!load_structure(input, structure) && mode != run_type::SYMM)
        return false;

    io::deserialize_basic(input,
        "add_hydrogen", add_hydrogen,
        "log_filename", log_filename,
        "n_threads", n_threads,
        "symm", symm_settings,
        "atac", atac_settings,
#ifdef SP2_ENABLE_LAMMPS
        "lammps", lammps_settings,
#endif // SP2_ENABLE_LAMMPS
#ifdef SP2_ENABLE_PHONOPY
        "phonopy", phonopy_settings,
#endif // SP2_ENABLE_PHONOPY
        "relax", relaxation_settings
    );

    return true;
}
