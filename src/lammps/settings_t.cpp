#include "lammps/settings_t.hpp"

using namespace sp2;

void lammps::lammps_settings_t::serialize(Json::Value &output) const
{
    io::serialize_basic(output,
        "n_tasks", n_tasks,
        "compute_lj", compute_lj,
        "compute_torsion", compute_torsion,
        "log_output", log_output,
        "sigma_scale", sigma_scale);
}

bool lammps::lammps_settings_t::deserialize(const Json::Value& input)
{
    return io::deserialize_basic(input,
        "n_tasks", n_tasks,
        "compute_lj", compute_lj,
        "compute_torsion", compute_torsion,
        "log_output", log_output,
        "sigma_scale", sigma_scale);
}
