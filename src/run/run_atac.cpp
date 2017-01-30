#include "run/run_types.hpp"

#include "common/util/mpi.hpp"
#include "common/minimize/minimize.hpp"

#include "atac/system_control_t.hpp"
#include "lammps/lammps_interface.hpp"
#include "airebo/system_control_t.hpp"

#include <memory>

using namespace std;
using namespace sp2;

int sp2::run_atac(const Json::Value &config, MPI_Comm comm)
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

    // check for/load minimization settings
    atac_settings_t atac_set;
    if (!io::deserialize_field(config, atac_set, "atac"))
        return EXIT_FAILURE;

    // mpi settings
    util::mpi_group_t mpi_info(0, comm);

    // pointer to hold our system once its constructed
    unique_ptr<atac::atac_subsys_t> sys(nullptr);


////////////////
// set up the system
////////////////

    // set up system depending on the input potential type
    switch (potential)
    {
    case potential_type::LAMMPS: {

        // load LAMMPS settings
        lammps_settings_t lmp_set;
        if (!io::deserialize_field(config, lmp_set, "lammps"))
            return EXIT_FAILURE;

        // setup mpi sub-communicators
        mpi_info = util::mpi_group_t(lmp_set.n_tasks, comm);

        // create the lammps system interface object
        unique_ptr<lammps::system_control_t> lmp_ptr(
            new lammps::system_control_t(mpi_info));

        // initialize with the input structure and other settings
        lmp_ptr->init(structure, lmp_set);

        // pass the pointer back out for minimization
        sys = move(lmp_ptr);
        break;
    }
    case potential_type::REBO: {
        // no settings to load for rebo, so just make it and pass it back out
        unique_ptr<rebo::system_control_t> rebo_ptr(
            new rebo::system_control_t());

        rebo_ptr->init(structure);
        sys = move(rebo_ptr);
        break;
    }
    default:
        return EXIT_FAILURE;
    }

////////////////
// create the controller object with the sub-system (which changes
// based on potential) and do the run
////////////////

    atac::system_control_t atac_sys(mpi_info);
    atac_sys.init(move(sys), atac_set);

    // do an initial minimization to avoid accepting
    // the first attempted mutation if the structure isn't
    // relaxed properly
    minimize_settings_t min_set;
    min_set.iteration_limit = 1000;
    min_set.output_level = 0;
    min_set.value_tolerance = 1e-4;

    minimize::acgsd(atac_sys, min_set);

    // output initial structure
    if (mpi_info.world_rank == 0)
        atac_sys.write_output("mutate.xyz");

    // iterate until the iteration limit is reached or an error is encountered
    while (atac_sys.iterate()) {}

    return EXIT_SUCCESS;
}
