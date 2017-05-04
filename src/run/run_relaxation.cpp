#include "run/run_types.hpp"
#include "common/minimize/minimize.hpp"

#include "lammps/lammps_interface.hpp"
#include "airebo/system_control_t.hpp"

int sp2::run_relaxation(const run_settings_t &config, MPI_Comm)
{
    if (config.relaxation_settings.output_file.empty())
    {
        std::cout << "No output file name given for structural relaxation."
                  << std::endl;

        return EXIT_FAILURE;
    }

    auto run_relax = [&](auto&& sys) {
        // minimize
        minimize::acgsd(sys.get_diff_fn(), sys.get_position(),
            config.relaxation_settings.minimize_settings);

        // write output
        if (!io::write_structure(config.relaxation_settings.output_file,
                sys.get_structure()))
        {
            return EXIT_FAILURE;
        }

        return EXIT_SUCCESS;
    };

    switch (config.potential)
    {
    case potential_type::LAMMPS: {
        lammps::system_control_t sys(config.structure, config.lammps_settings);

        return run_relax(sys);
    }
    case potential_type::REBO: {
        airebo::system_control_t sys(config.structure);

        return run_relax(sys);
    }
    default:
        return EXIT_FAILURE;
    }
}
