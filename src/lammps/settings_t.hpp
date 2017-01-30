#ifndef SP2_LAMMPS_SETTINGS_T_HPP
#define SP2_LAMMPS_SETTINGS_T_HPP

#include "common/json/json.hpp"

namespace sp2 {

struct lammps_settings_t : public io::json_serializable_t
{
    int n_tasks;            ///< number of MPI tasks per instance of LAMMPS (0 for all available)
    bool compute_lj,        ///< flag to compute the lennard-jones forces (van der Waals)
        compute_torsion,   ///< flag to compute torsion for C-C bonds
        log_output;        ///< whether lammps should write to stdout or not
    double sigma_scale;     ///< range of lj forces is this times 3.4 angstroms

    lammps_settings_t();

    bool serialize(Json::Value& output) const;
    bool deserialize(const Json::Value& input);
};

} // namespace sp2

#endif // SP2_LAMMPS_SETTINGS_T_HPP
