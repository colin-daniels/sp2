#ifndef SP2_LAMMPS_SETTINGS_T_HPP
#define SP2_LAMMPS_SETTINGS_T_HPP

#include "common/json/json.hpp"

namespace sp2 {
namespace lammps {

struct lammps_settings_t : public io::json_serializable_t
{
    int n_tasks = 0;             ///< number of MPI tasks per instance of LAMMPS (0 for all available)
    bool compute_lj = true,      ///< flag to compute the lennard-jones forces (van der Waals)
        compute_torsion = false, ///< flag to compute torsion for C-C bonds
        log_output = false;      ///< whether lammps should write to stdout or not
    double sigma_scale = 3.0;    ///< range of lj forces is this times 3.4 angstroms

    void serialize(Json::Value &output) const;

    bool deserialize(const Json::Value &input);
};

} // namespace lammps
} // namespace sp2

#endif // SP2_LAMMPS_SETTINGS_T_HPP
