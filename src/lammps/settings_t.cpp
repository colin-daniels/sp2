#include "lammps/settings_t.hpp"

using namespace sp2;

lammps_settings_t::lammps_settings_t() :
    n_tasks(0),
    compute_lj(true),
    compute_torsion(false),
    log_output(false),
    sigma_scale(3.0) {}

bool lammps_settings_t::serialize(Json::Value& output) const
{
    io::serialize_basic(output,
        "n_tasks", n_tasks,
        "compute_lj", compute_lj,
        "compute_torsion", compute_torsion,
        "log_output", log_output,
        "sigma_scale", sigma_scale);
    return true;
}

bool lammps_settings_t::deserialize(const Json::Value& input)
{
    io::deserialize_basic(input,
        "n_tasks", n_tasks,
        "compute_lj", compute_lj,
        "compute_torsion", compute_torsion,
        "log_output", log_output,
        "sigma_scale", sigma_scale);
    return true;
}
