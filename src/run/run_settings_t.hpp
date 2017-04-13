//
// Created by cc on 11/25/16.
//

#ifndef SP2_RUN_SETTINGS_T_HPP
#define SP2_RUN_SETTINGS_T_HPP

#include "atac/settings_t.hpp"
#include "relaxation/relaxation_settings.hpp"
#ifdef SP2_ENABLE_PHONOPY
#include "phonopy/phonopy_settings.hpp"
#endif // SP2_ENABLE_PHONOPY

#ifdef SP2_ENABLE_LAMMPS
#include "lammps/settings_t.hpp"
#endif // SP2_ENABLE_LAMMPS

#include "symm/symm_settings_t.hpp"
#include "common/io/structure.hpp"
#include "common/structure_t.hpp"
#include "common/enums.hpp"
#include "common/json/json_serializable_t.hpp"

namespace sp2 {

struct run_settings_t : public io::json_serializable_t
{
    run_type mode = run_type::NONE;
    potential_type potential = potential_type::LAMMPS;

    structure_t structure;

    bool add_hydrogen = false;
    /// number of threads to use (when applicable):
    ///     0 for off
    ///     -1 for std::thread::hardware_concurrency
    int n_threads = 0;
    std::string log_filename = "progress.log";

    symm::symm_settings_t symm_settings;
    atac::atac_settings_t atac_settings;

#ifdef SP2_ENABLE_LAMMPS
    lammps::lammps_settings_t lammps_settings;
#endif // SP2_ENABLE_LAMMPS
#ifdef SP2_ENABLE_PHONOPY
    phonopy::phonopy_settings_t phonopy_settings;
#endif // SP2_ENABLE_PHONOPY

    relaxation_settings_t relaxation_settings;

    bool serialize(Json::Value &output) const;
    bool deserialize(const Json::Value &input);
};

} // namespace sp2


#endif // SP2_RUN_SETTINGS_T_HPP
