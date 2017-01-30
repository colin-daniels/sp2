#include "run/run_types.hpp"

#include "common/util/mpi.hpp"
#include "common/minimize/minimize.hpp"

#include "lammps/lammps_interface.hpp"
#include "airebo/system_control_t.hpp"

using namespace std;
using namespace sp2;

int sp2::run_minimize(const Json::Value &config, MPI_Comm comm)
{
////////////////
// read potential type (e.g. LAMMPS, REBO) and structure
////////////////

    potential_type potential;
    if (!load_potential_type(config, potential))
        return EXIT_FAILURE;

    structure_t structure;
    if (!load_structure(config, structure))
        return EXIT_FAILURE;

////////////////
// read primary settings for the run from the config
////////////////

    minimize_settings_t min_set;
    if (!io::deserialize_field(config, min_set, "minimize"))
        return EXIT_FAILURE;

    // mpi settings, note that input numbers of tasks are ignored for all
    // potential types, since there's no reason to use mpi groups
    util::mpi_group_t mpi_info(0, comm);

    // pointer to hold our system once its constructed
    unique_ptr<minimize::differentiable_t> sys(nullptr);

////////////////
// set up the system depending on the input potential type
////////////////

    switch (potential)
    {
    case potential_type::LAMMPS: {

        // load LAMMPS settings
        lammps_settings_t lmp_set;
        if (!io::deserialize_field(config, lmp_set, "lammps"))
            return EXIT_FAILURE;

        if (mpi_info.world_size > lmp_set.n_tasks &&
            lmp_set.n_tasks != 0)
            cout << "Warning, task grouping for LAMMPS is ignored when "
                 << "running minimization only (all processors are "
                 << "automatically used)." << endl;

        // create the lammps system interface object, note that minimization
        // ignores number of tasks input in settings and just uses all the
        // available tasks
        mpi_info = util::mpi_group_t(numeric_limits<int>::max(), comm);
        unique_ptr<lammps::system_control_t> lmp_ptr(
            new lammps::system_control_t(mpi_info));

        // initialize with the input structure and other settings
        lmp_ptr->init(structure, lmp_set);

        // pass the pointer back out for minimization
        sys = move(lmp_ptr);
        break;
    }
    case potential_type::REBO: {
        if (mpi_info.world_size > 1)
        {
            if (mpi_info.world_rank != 0)
                return EXIT_SUCCESS;

            cout << "Warning, the REBO potential mode currently only uses one "
                 << "thread, multiple tasks do nothing in minimization mode."
                 << endl;
        }

        // no settings to load for rebo, so just make it
        unique_ptr<rebo::system_control_t> rebo_ptr(
            new rebo::system_control_t());

        rebo_ptr->init(structure);

        // and pass it back out
        sys = move(rebo_ptr);
        break;
    }
    default:
        return EXIT_FAILURE;
    }

////////////////
// perform the minimization
////////////////

    // run the minimization routine
    minimize::acgsd(*sys, min_set);

    if (!min_set.output_filename.empty() && mpi_info.world_rank == 0)
        sys->write_output(min_set.output_filename);

    return EXIT_SUCCESS;
}
